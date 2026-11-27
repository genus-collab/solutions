# Frame-B formula and refined conifold gap for local P2 x C^2 on the fibrewise torus: proof of record

File: hodge_vertex/lp2_explore/FRAMEB_PROOF.md. Written 2026-11-02 (round 6c); edited in place in rounds 6c, 6d, 6e and 6f. Round 6e also writes H^* for A^* in the injectivity of p^* (A1.7, A3.3 Colour isolation); the projection-formula argument is unchanged.

**Provenance.** Assembled from the prover reports FB_A, FB_B, FB_C and the referee reports refFB_A, refFB_A2, refFB_B, refFB_C, refFB_C2 in hodge_vertex/lp2_explore/round6_scripts/round6b_result.json. All five referee verdicts are "sound with minor fixes", with no major item. Each referee fix is applied here or explicitly declined (Appendix, s.A.5). The fixes that change mathematics were re-derived by hand before being applied. The step numbering A1.1, A2.4, B4, C Step 6, ... follows the reports. **The referee passes of this assembled text are recorded at the marker ROUND6C-REFEREES in s.0.4 of the draft. They are AI passes; no human referee has read this proof.** The JSON keys A, A2, B, C, C2 are refFB_A, refFB_A2, refFB_B, refFB_C, refFB_C2. Minor items of the round-6c referee reports are applied in place where they could be re-derived; A2.8 is the one place where a proof route changed. Minor items of the round-6d reports ref28hand, ref28num, refNew, refCon2, lit, coldF, coldB are applied in place after checking against this text. The one major item of round 6d (report lit, citation level) is item O1 below. Its citation was repaired in round 6e, after the verification report verify:ccit, whose corrections to lit's text are applied. Items F1-F8 of the round-6f report refQDE (round6_scripts/round6f_result.json), the second reader of the QDE derivation of s.1.2, are applied in s.0.2, s.1.2 and s.5.1; no later pass has read the inserted text.

**Tags.** **[P]** proved here, to all orders. **[C]** cited; s.1 says what was opened. **[E]** finite-order evidence; never a proof step. "P mod X" means proved here given X, and nothing beyond X.

**Rulings respected.** BS genus-0 mirror and oscillating-integral statements on T_F may be used. BS Thm cfasym is acceptable as a source of F2 and is not used for it here (F2 follows from the closed form, draft Cor 2.9). BS/BFGW finite generation and algebraic independence are acceptable. Topological recursion on the local-P2 mirror curve is not used. BS rmk:refmm, BKMP/FLZ remodelling, HK10, Coates-Iritani 2018, BFGW Conj 1.9, BS prop:cfgraphsum, thm:cfasym and lem:cfpole are not inputs. No finite-order computation is a proof step. Scope is the fibrewise torus only, per colouring. Stage A uses eq:J=I and the genus-0 QDE at generic alpha on T_max, not on T_F. eq:J=I at T_max is CCIT15 Theorem 31 (s.1.2), a published all-orders genus-0 theorem whose proof uses localisation and the Lagrangian cone, with no Eynard-Orantin TR, remodelling or other forbidden input. BS's own citation MR2486673 covers T_F only. The genus-0 QDE at T_max is derived in s.1.2. Ruling Q1 names statements on T_F. This use of CCIT15 Theorem 31 is outside the literal wording of Q1; the user confirmed it on 2026-11-02.

**Item O1 (citation repaired in round 6d-6e; citation level; no step of Stages A-C changes; s.1.2 gains a derivation of the genus-0 QDE at T_max, which was previously taken as a cited whole).** The chain needs eq:J=I at generic alpha on T_max. *What BS cites.* BS 2658 cites MR2486673 = Coates, Commun. Math. Phys. 287 (2009) 1071-1108, arXiv:0810.2200v1. *What it covers.* Its torus covers the trivial action on P2, so that its Prop 4.1, Cor 4.2 are eq:J=I on T_F only (lambda = E, p = H + E/3). Stage A cannot be moved to T_F, because the Q = 0 product is not semisimple there. *What is now cited.* The source of record is CCIT15 Theorem 31 (s.1.2), a published genus-0 theorem whose proof uses localisation and the Lagrangian cone, with no Eynard-Orantin TR, remodelling or other forbidden input. *How it was verified.* Report lit (round 6d) found the defect and proposed the source. An independent reader (report verify:ccit, round 6e, an AI pass) re-read the statement and the definitions, did the substitution, the mirror-map step and the comparison with A1.2 by hand, and ran an exact finite-order check, which is evidence only (s.A.1). The proof of Theorem 31 was read for structure and not re-derived, and the journal print was not read. *QDE source.* No numbered theorem at T_max was found. The genus-0 QDE at T_max is derived in s.1.2 from Pandharipande's Lemma 1 and A1.4. The T_max use is outside the literal wording of ruling Q1, which names statements on T_F. The user confirmed it on 2026-11-02, accepting CCIT15 Theorem 31 at generic weights under that ruling. The statements below that rest on eq:J=I are: A2.1, A2.4-A2.6, A2.9, A3, parts (A2), (A3) and the row form of the graph sum of Theorem FB_A; Theorem FB_B except (b); Proposition 13.2_P2; Corollaries G and G_{k,h}; the discharge paragraph at the end of s.2. Not affected: A1, A2.2, A2.3, A2.7, A2.8, A3.2, the first paragraph of A3.3 (the part A1.7 uses), A3.4, the comparison statement of A2.9 for two flat solutions (only its application to the R~ of A2.4 rests on eq:J=I), the identity Tbar_g[𝖥] = G_g[𝖥] of Theorem FB_C (C1), (C2), Lemma (t_cf), and draft Theorems LC_P2, R_P2, Corollary SV_P2, Theorem VI as statements about the decoration sum; also draft SP1-SP5, Lemmas 2.1, 2.3-2.5, Theorem A, Cor 2.2, Cor 2.8, Props 2.10, 2.12. This partition is a dependency reading of the current text (round-6d fix step), re-read once by the round-6e audit, an AI read following the citations in A1-A3, B5, C Steps 0-2 and C Step 6. No referee report has checked it.

Contents. 0 Summary. 1 Inputs. 2 Stage A: the CohFT, the R-matrix and its constants at generic weights. 3 Stage B: specialisation to T_F. 4 Stage C: frame change, Cf, t_cf, the gap. 5 Errata for BS and the draft. Appendix: evidence, and the disposition of the referee items.

---

## 0. Summary

### 0.1 Notation

The draft is hodge_vertex/lp2_explore/LP2_REGULARITY_PROOF_draft.md; "s.n of the draft" refers to it. BS is BS_v2_2026-11-02/main.tex, cited by label and line.

**Weights and fields.**
- E := eps1 + eps2. BS's function is always written E(y), with its argument. The Hodge bundle is written 𝔼. Teleman's operator is written 𝐄.
- wt_i := {3 alpha_i, alpha_j - alpha_i, alpha_k - alpha_i, -eps1, -eps2}, with W_i^0 := the product of wt_i. T_F is alpha* = (E/3)(1,1,1), and w := (E,-eps1,-eps2).
- Kw := Q(alpha_0,alpha_1,alpha_2,eps1,eps2) and K~w := Kw(sqrt W_0^0, sqrt W_1^0, sqrt W_2^0).
- K := Q(eps1,eps2), K_y := K(y), L := K(t) with t^3 = y, and L_omega := K(omega)(t) with omega = e^{2 pi i/3}. (The report FB_B writes K, L for K_y, L. The report FB_C writes L for L_omega.)
- For a multiset v: h_v(z) := exp(-Sum_k [B_{2k}/(2k(2k-1))] p_{2k-1}(1/v) z^{2k-1}) = prod_{a in v} Theta(-z/a), where Theta(x) := exp Sum_{m>=1} B_{2m} x^{2m-1}/(2m(2m-1)). Further Hodge_v := prod_a Lambda^vee(v_a)/v_a^g, with Lambda^vee(u) := Sum_i (-1)^i lambda_i u^{g-i}. BS's c_j satisfies e^{-c_j} = h_{wt_j}.
- h_w(z) = Theta(-z/E) Theta(z/eps1) Theta(z/eps2), and Omega_w := Hodge_w.

**Base.**
- theta := theta_y = y d/dy, Delta_c := 1 + 27y = xi, X := 1/xi.
- The mirror map is Q = y exp(3 Sum_{d>=1} (3d-1)!/(d!)^3 (-y)^d). It is a unit times y with coefficients in Q and independent of alpha, so Kw[[Q]] = Kw[[y]]. A := theta_y log Q = Sum_d a_d (-y)^d with a_d = (3d)!/(d!)^3, theta_Q = A^{-1} theta_y, F := theta log A.
- Lhat(Y) := prod_j (Y + alpha_j) + 27 y Y^3 and W(Y) := e_1(alpha) Y^2 + 2 e_2(alpha) Y + 3 e_3(alpha).
- 𝖥 and 𝖠 denote indeterminates; F, A are their values.

**Rings of Stage B.**
- 𝔸 := Q[alpha,eps,y], with the prime 𝔮 := (alpha_i - E/3, i = 0,1,2). O := 𝔸_𝔮, O/𝔮O = K_y. O_w := Q[alpha,eps] localised at 𝔮 ∩ Q[alpha,eps].
- 𝒞 := O[Y]/(Lhat). Reduction rho: O_w[[y]] -> K[[y]] is alpha_i -> E/3 coefficientwise.

**Rows.**
- Rows at generic alpha: R~, R, P_a^{(i)} := Delta_i^{1/2} R_{ai}. Algebraic row: Phi, P_0, 𝐓, 𝐐 of Theorem FB_B(a).
- Rows on T_F: t_i := omega^i t, Y_i := -(E/3)/(1+3t_i), v_i := 3/(E+3Y_i) = (1+3t_i)/(E t_i), Delta_i := 3 eps1 eps2 E/v_i^2; S(v,z) := Sum_{l>=0} r_l (vz)^l/prod_{j=1}^{l+1}(1+jz/E), with r_l from (Rec) of s.1.3 of the draft.

**Graph sums.** For g >= 2, G_g is the decoration sum (0.1) of the draft. 𝒯bar_g(F) is defined in Theorem FB_B(d). Tbar_g[𝖥] and G_g[𝖥] are defined in the Notation of s.4, and Tbar_g[𝖥] is the reduction of the polynomial 𝒯_g of Theorem FB_B(c) (this identification is B6 and rests on eq:J=I; the identity of Theorem FB_C (C1) does not use it).

### 0.2 Statements and status

**Theorem FB_A (generic weights).** [P mod eq:J=I (CCIT15 Thm 31, s.1.2; item O1) and the genus-0 QDE (derived in s.1.2; eq:SHb is derived from them in A2.1), on which (A2), (A3) and the row form of the graph sum rest; Teleman Theorem (2), Prop 8.7(ii), and 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss); PPZ Defs 2.2, 2.5, 2.13, Prop 2.14; Graber-Pandharipande; the fundamental-class and splitting axioms of the virtual class, with Coates-Givental App. 1 (iv)-(vi) for the twist (Behrend-Manin not opened); Mumford's theorem as in Faber-Pandharipande s.1.1; Faber-Pandharipande Theorem 4. Part A3.5 additionally mod draft SP1-SP5 and Lemmas 2.1, 2.3-2.5. Two referee passes on the report form; passes on the assembled text as recorded at the marker ROUND6C-REFEREES in s.0.4 of the draft. A2.8 was rewritten in round 6c with a direct derivation in place of the report's route through Teleman's Prop 7.7; that derivation has had two independent hand re-derivations, recorded in A2.8.]

(A1) Let Omega_{g,n} := Sum_d Q^d st_*([Mbar_{g,n}(P2,d)]^vir_T cap prod ev^*gamma_k cap e_T(R pi_* f^* N)^{-1}), with N = O(-3)_{3 alpha} + O_{-eps1} + O_{-eps2}. This equals BS's localisation definition. It is a CohFT with flat unit over Kw[[Q]], with pairing eta(a,b) = Sum_i a|_{P_i} b|_{P_i}/W_i^0. It satisfies the divisor equation p_* Omega_{g,n+1}(gamma, H) = theta_Q Omega_{g,n}(gamma). It is semisimple at Q = 0, with idempotents e_i(Q) in V (x) Kw[[Q]] and Delta_i(Q = 0) = W_i^0. Each GW^T_{g,d} (d >= 1, or d = 0 and g >= 2) lies in O_w. So it is regular on the diagonal, and its value at alpha* is the T_F integral Int_{[Mbar(P2,d)]^vir} e_{T_F}(R pi_* f^* N)^{-1}. For any embedding of K~w((Q)) in ℂ, Teleman gives Omega = R.omega (PPZ Def 2.13), with omega(e_i^n) = Delta_i^{g-1}, for a unique symplectic R, and R has coefficients in K~w((Q)).

(A2) Put s_{aj} := (W_j^0)^{-1/2} S(H^a)|_{P_j}. Define u_j by theta_Q u_j = lambda_j, where lambda_j is the eigenvalue of H* with lambda_j(0) = -alpha_j, normalised by u_j + alpha_j log y -> 0 at y = 0. Then theta_y u_j = Y_j, where Y_j in Kw[[y]] is the root of Lhat with Y_j(0) = -alpha_j. Moreover R~_{aj} := e^{-u_j/z} s_{aj} lies in K~w[[y]][[z]], R~|_{z=0} = Psi (Psi_{aj} = f_j^*(H^a)), R~|_{y=0} = Psi|_{y=0} (no z-dependence), Delta_i = 3 eps1 eps2 W(Y_i), and Sum R~_{ai}(z) eta^{ab} R~_{bj}(-z) = delta_ij. Further R~_{1j} = (z theta_y + Y_j) R~_{0j}/A and R~_{2j} = (z theta_y + Y_j - e_1 E(y)) R~_{1j}/D (eq:RHb). With the dictionary eta(H^a, R(z) f_j) = R_{aj}(z), Teleman's R satisfies R_{aj}(z) = R~_{aj}(z) exp(Sum_k a_k^{(j)} z^{2k-1}) with constants a_k^{(j)} in K~w.

(A3) BS (Normalisation) is a theorem: a_k^{(j)} = -[B_{2k}/(2k(2k-1))] p_{2k-1}(1/wt_j). That is, R_{aj} = R~_{aj} e^{-c_j} (eq:canR), R_{aj}|_{Q=0} = Psi_{aj}|_{Q=0} h_{wt_j}(z), and in End(V): R|_{Q=0} = diag(h_{wt_j}(z)) on the f_j. For n >= 1, Omega(e_{i_1},...,e_{i_n})|_{Q=0} = delta_{i_1=...=i_n=i} (W_i^0)^{g-1} Hodge_{wt_i}. CCIT Prop 6.9 is not needed: with Asym as in BS (normalised by sqrt det Hess(-W)), Asym_{q_i^cr}[e^{W/z} dOmega] = c N R_{0i}, with N = sqrt(eps1 eps2) Theta(-z/eps1) Theta(-z/eps2), c^2 = -1, and c universal up to the branch of the square root. With (det Hess W)^{-1/2} in place of (det Hess(-W))^{-1/2}, c = 1. Structurally, Omega^Z_{g,n} = (eps1 eps2)^{g-1} Hodge_{(-eps1,-eps2)} Omega^{K_P2}_{g,n}, so R^Z = R^{K_P2} h_{(-eps1,-eps2)}, and Asym = c R^{K_P2}_{0i}, with R^{K_P2} taken in its own normalised frame, so that c = 1 in the (det Hess W)^{-1/2} convention. (R^Z_{0i} = (eps1 eps2)^{-1/2} R^{K_P2}_{0i} h_{(-eps1,-eps2)} and Theta(-z/eps) Theta(z/eps) = 1.)

(Graph sum, g >= 2.) In Kw[[Q]] (each coloured term lies in Kw[[Q]], with no square roots and no cube roots at generic alpha; the identity is for the total),

    F^T_g := Sum_{d>=0} Q^d GW^T_{g,d}
           = Sum_Gamma |Aut Gamma|^{-1} Sum_{c: V -> {0,1,2}} Sum_dec prod_v [Delta_{c_v}^{g_v-1} (1/m_v!) Int_{Mbar_{g_v,n_v+m_v}} prod psi^{a_h} prod psi^{d_k} prod_k Dil_{c_v}(d_k)] prod_e E_{cc'}(a_e,b_e),

where P_a^{(i)} := Delta_i^{1/2} R_{ai}, Dil_i(d) := (-1)^d [z^{d-1}] P_0^{(i)} (d >= 2), and E_ij(psi,w) := [delta_ij Delta_i - Sum P_a^{(i)}(-psi) eta^{ab} P_b^{(j)}(-w)]/(psi + w). Equivalently, F^T_g is the R~-sum with vertex classes Delta_{c_v}^{g_v-1} Hodge_{wt_{c_v}}.

F^T_g includes d = 0. BS's GW_g as printed (line 1277) has d > 0 only, as do the correlators of eq:gwprim (line 2446).

**Theorem FB_B (specialisation to T_F; replaces draft Prop 2.11).** [P mod Theorem FB_A (A1.2, A2, A3.5, graph sum); draft SP1-SP5, Lemmas 2.1, 2.3-2.5, Theorem A, Cor 2.2, Prop 2.10, Prop 2.12; Graber-Pandharipande; c(𝔼)c(𝔼^vee) = 1 (Mumford (5.4), (5.5); FP s.2 before (22)), for the constant term. One referee pass in report form, conditional on Stage A; later passes as recorded at the marker ROUND6C-REFEREES. Through Theorem FB_A, also P mod eq:J=I (CCIT15 Thm 31, s.1.2), except (b).]

Let F_g^{T_F}(y) := Sum_{d>=0} Q(y)^d GW^{T_F}_{g,d} in K[[y]], g >= 2.

*Algebraic row.* Put Phi(alpha,Y;z) := exp((z/2) G^{ab} d_a d_b) exp(z^{-1} Sum_{m>=3} T_m(xi))|_{xi=0}, with G = -H^{-1}, H = diag(alpha_j+Y) - (Y/3) J (J all ones), and T_m = Sum_j (alpha_j+Y) xi_j^m/m! - 3Y (Sum_j xi_j/3)^m/m!. Then det H = W/3, Phi_0 = 1, and [z^n] Phi lies in Q[alpha,Y,W^{-1}]. Put P_0 := Theta(z/eps1) Theta(z/eps2) Phi, 𝐓 := z theta + Y - (z/2) theta log W, 𝐐 := 𝐓 P_0.

*Rows on T_F (T_F limit of the canonical rows).* Pbar_0^{(i)} := h_w(z) S(v_i,z), Qbar := 𝐓bar Pbar_0 with 𝐓bar := z theta + Y_i - (z/2) theta log Delta_i, Pbar_1 := Qbar/A, Pbar_2 := Delta_c A [𝐓bar Qbar - z F Qbar - E (X/2)(A^{-2} - 1) Qbar]. These are the large-radius rows of draft Def 1.3 at c1 = A, multiplied by Delta^{1/2}. Edges and dilatons are as in Stage A with etabar^{ab} = eta^{ab}|_{alpha_i=E/3}.

(a) *Regularity of factors, all orders in z.* 𝒞 is free of rank 3 over O and local, with maximal ideal 𝔮𝒞 and residue field L, via pi: Y -> -(E/3)/(1+3t), y -> t^3. The elements Y, W, Y+alpha_j, Lhat'(Y) = -W/Y are units of 𝒞. theta extends uniquely to 𝒞, with theta Y = -Y prod_j(Y+alpha_j)/W, reducing to (t/3) d/dt. P_0, 𝐐, 𝐓𝐐 have z-coefficients in 𝒞: no square root and no discriminant in the denominators, only powers of W, eps1, eps2. For every colour j at generic alpha, Delta_j^{1/2} R_{0j} = iota_j(P_0), where iota_j: 𝒞 -> Kw((y)), Y -> Y_j(y). Hence Delta_j, Delta_j^{-1} and Dil_j(d) are iota_j of elements of 𝒞, and, with 𝖥 an indeterminate,

    (Pair) [z^M w^b]<P(z),P(w)> = eta^{00} P_{0,M} (x) P_{0,b} - 3 e_1 eps1 eps2 𝐐_M (x) 𝐐_b
             - 3 eps1 eps2 Delta_c (𝐐_M (x) (𝐓𝐐)_b + (𝐓𝐐)_M (x) 𝐐_b)
             + 3 eps1 eps2 Delta_c 𝖥 (𝐐_{M-1} (x) 𝐐_b + 𝐐_M (x) 𝐐_{b-1}),

with eta^{00} = 3 e_3 eps1 eps2. A has cancelled, and E_ij(a,b) = -Sum_{k=0}^b (-1)^k (-1)^{a+b+1} [z^{a+k+1} w^{b-k}]<P(z),P(w)> lies in (iota_i (x) iota_j)(𝒞 (x)_O 𝒞) + F (iota_i (x) iota_j)(𝒞 (x)_O 𝒞), linear in F. Equivalently, d_𝖥 E_ij(psi,w) = 3 eps1 eps2 Delta_c 𝐐^{(i)}(-psi) 𝐐^{(j)}(-w). The symplectic identity and the delta form of the edge hold with 𝖥 free (B3). pi carries all of this to the T_F data: pi(Phi) = Theta(-z/E) S(v,z), pi(P_0) = h_w S, at each z-order.

(b) *GW side (BS line ~4056, proved in a stronger form).* For d >= 1 and g >= 0, eps1 eps2 GW^T_{g,d} lies in Q[alpha,eps], is homogeneous of degree 2g, and GW^T_{g,d} = Int_{[Mbar_g(P2,d)]^vir_T} e_T(R^1 pi_* f^* O(-3)_{3 alpha}) Lambda^vee(-eps1) Lambda^vee(-eps2)/(eps1 eps2). For d = 0 and g >= 2, 3 e_3 eps1 eps2 GW^T_{g,0} lies in Q[alpha,eps], using (-3H)^{-1} = (H^2 + e_1 H + e_2)/(3 e_3) in H_T(P2)[e_3^{-1}], with H the equivariant class, H|P_i = -alpha_i. Both lie in O_w, and their image under alpha_i -> E/3 is GW^{T_F}_{g,d} = Int_{[Mbar_g(P2,d)]^vir} e_{T_F}(R pi_* f^* N)^{-1}, the T_F localisation definition. Constant term, where H' is the non-equivariant hyperplane class of P2 and Int_{P2} H'^2 = 1:

    GW^{T_F}_{g,0} = Int_{Mbar_g} [9 lambda_g lambda_{g-2} P(0) - 3 lambda_g lambda_{g-1} P'(0)],
    P(H') := Lambda^vee(E-3H') Lambda^vee(-eps1) Lambda^vee(-eps2)/((E-3H') eps1 eps2).

For g = 2 this is E^2/640 [P, by hand]. For g = 3 this is 9 E^4 Int lambda_1 lambda_2 lambda_3 = E^4/161280 [P, by hand, with FP Theorem 4; written in B5]. No closed form is claimed for g >= 4.

(c) *Generic alpha, algebraic form.* For each stable Gamma of genus g >= 2 without legs and each decoration, the term of Stage A's graph sum at colouring c is (prod_v iota_{c_v}) applied to an element T_{Gamma,dec} of 𝒞^{(x)_O |V|}[𝖥]. So 𝒯_g := Sum_Gamma |Aut|^{-1} Sum_dec (Tr_{𝒞/O})^{(x)|V|}(T_{Gamma,dec}) lies in O[𝖥], deg_𝖥 <= 3g-3, and Stage A's identity reads F^T_g = 𝒯_g(F) in Kw[[y]], where the sum over colourings is the multi-trace.

(d) **The specialisation identity.** In K((y)),

    F_g^{T_F}(y) = Sum_Gamma |Aut Gamma|^{-1} Sum_{c: V -> {0,1,2}} Sum_dec prod_v [Delta_{c_v}^{g_v-1} (1/m_v!) Int psi^a psi^d prod_k Dilbar_{c_v}(d_k)] prod_e Ebar_{c c'}(a_e,b_e) =: 𝒯bar_g(F),

the point-frame graph sum for the rows Pbar. Each coloured term lies in K(omega)[t^{+-1},(1+3 omega^i t)^{-1} : i][F] and has poles at t = 0 in general. The colour sum lies in K[y^{+-1},Delta_c^{-1}][F]. The left side is a formal series in K[[y]]; its convergence near y = 0 is a consequence of (d), because the right side is algebraic in F and y. Continuation to the universal cover of ℂ \ {0,-1/27} means continuation of the right side, and it is not a statement about GW series at large y. The T_F CohFT is not semisimple at Q = 0, and Teleman is never applied on T_F.

(e) *Consequences.* F_g^{T_F} lies in K[y^{+-1},Delta_c^{-1}][F] with deg_F <= 3g-3. No BS prop:fingen* is used, and neither is A, since A cancels in leg-free graphs. F is transcendental over ℂ(y) by an elementary monodromy argument, so the polynomial in 𝖥 is unique, and Cf(F_g^{T_F}) := 𝒯bar_g(F_cf) is well defined. Colouring by colouring, the same holds over ℂ(t).

(f) *Draft Prop 2.11.* The family eq:adelta and BS rational tori are unnecessary. The step "if delta divides Q some phi_{k,m} blows up" is deleted, not repaired: non-divisibility is obtained a priori as den not in 𝔮 from (a). Algebraic independence is not used for (d).

**Theorem FB_C (C1): frame change, for all g >= 2.** [P mod draft Prop 2.10, A3.2 (Mumford as in FP), PPZ Defs 2.2, 2.5, 2.13, Prop 2.14. Stages A, B are needed only to identify Tbar_g[𝖥] with the reduction of 𝒯_g, B6, not for this identity.] In L_omega[𝖥], Tbar_g[𝖥] = G_g[𝖥]. This is a polynomial identity in 𝖥. It holds for the total, not graph by graph or colouring by colouring. More precisely, Omega_w = h_w.1 for h_w = Theta(-z/E) Theta(z/eps1) Theta(z/eps2), with PPZ Def 2.13 conventions; the rows of R' h_w are h_w(z) times the rows of R', so R' = R_{T_F} h_w^{-1}. The exponent -1 is derived. h_w is the same at all colours, does not depend on t, and h_w(z) h_w(-z) = 1.

**Theorem FB_C (C2): dictionary (notation of s.4).** [P mod PPZ Defs 2.2, 2.5, 2.13 and draft Prop 2.10; at c1 = A also draft Prop 2.12.] T'(z) = z(1 - R'^{-1}(z)1) = Sum_i e_i z(1 - S(v_i,-z)) = Sum_i e_i Sum_{d>=2} Dil^_i(d) z^d, with Dil^_i(d) := (-1)^d [z^{d-1}] S(v_i,z) as in draft s.1.4, and with vanishing z^0 and z^1 terms. The code's q_t^{(i)} = -psi S(v_i,-psi) equals T'_i - psi, so -psi is exactly the dilaton shift (Dil^_i(1) = -1). The edge is (1.1) of the draft, and Ehat_ij[𝖥](psi,w) = Ehat_ij[0](psi,w) + 3 eps1 eps2 Delta_c 𝖥 Q_i(-psi) Q_j(-w), with Ehat as in (1.1) of the draft; in the point frame the 𝖥-term of Ebar carries the additional factor h_w(-psi) h_w(-w). Neither c1 nor A occurs in a leg-free term. Each decorated coloured term lies in K(omega)[t^{+-1},(1+3 omega^i t)^{-1}: i][𝖥]. At c1 = A, 𝖥 = F these are eq:RHb, eq:leafterms, eq:edgeterm (draft Prop 2.12).

**Definition (Cf), the reading of BS adopted.** BS 5148-5160 define Cf on generators by Cf(A) = A_cf, Cf(F) = F_cf. The lines fixing X and the Y_i (BS 5159, 5174) are commented out, and BS 5163-5167 contradicts them (s.5). The reading adopted here as the definition is: Cf is the identity on coefficients, Cf(X) = X, Cf(Y_i) = Y_i. It is forced by BS 5172 (D_cf = X/A_cf^2) and BS 5340 (the formula for Y_0^cf). Proposition 13.2_P2 is a theorem about this Cf. Formally:
- L~ := ℂ(eps1,eps2,eps3)(t) with eps3^2 = E, and t a formal variable; iota: G_{T_F} -> L~, Y_i -> Y_i, W_i^{1/2} -> eps3 (3Y_i+E)/3; 𝐆 := L~[𝖠^{+-1},𝖥]. The statement about the rows Cf(R_{ai}) needs sqrt(3 eps1 eps2) adjoined to L~. Leg-free terms do not.
- The target of ev_A is M(U)(eps1,eps2)(eps3), with M(U) the meromorphic functions on a small open U and a chosen branch of t on U. ev_A: 𝖠 -> A, 𝖥 -> F.
- ev_cf: 𝖠 -> A_cf, 𝖥 -> F_cf, t -> -(1-xi)^{1/3}/3, with target the same field built on meromorphic germs at xi = 0.
- Coloured Cf is ev_cf on the lift in 𝐆, with t a formal variable. On functions, Cf := ev_cf o ev_A^{-1}. For a single coloured function this depends on the branch of t on U (t -> omega t relabels colours). Only the colour-summed total is branch-independent.
- The value of Cf on GW: Cf(F_g^{T_F}) := ev_cf of the unique polynomial of Theorem FB_B(e).

**Proposition 13.2_P2 (d >= 0 reading).** [P mod Stages A, B; draft Props 2.10, 2.12; PPZ Defs 2.2, 2.5, 2.13, Prop 2.14 (C Steps 1, 2); and, for the clause F_cf = XX only, the logarithmic form in Lemma (t_cf) (DLMF 15.6.7, Gauss), not Barnes's second lemma. Through Stages A, B, also P mod eq:J=I (CCIT15 Thm 31, s.1.2).] For g >= 2, Cf(F_g^{T_F}) = G_g[F_cf], and this is G_g of (0.1) of the draft at c1 = C1cf. With BS's printed d > 0 convention, Cf(GW_g) = G_g[F_cf] - GW^{T_F}_{g,0}, with GW^{T_F}_{g,0} the constant in K of Theorem FB_B(b). BS's eq:F2ref includes the d = 0 term.

**Lemma (t_cf).** [P mod DLMF 15.6.7, Barnes's second lemma, Gauss's values of psi at 1/3, 2/3. Barnes's lemma confirmed only through a secondary statement, s.1.2.] On y in (-1/27,0) with log y = log(-y) + i pi: BS's t_cf = (sqrt3/(2pi))(l^2/2 + g l + g2 - pi^2/2), with l = log(-y), g = I_0^{[1]} - log y, g2 the series part of I_1^{[2]}; A_cf = -2F1(1/3,2/3;1;xi) = C1cf, F_cf = theta log A_cf = XX|_{C1cf}, and t_cf = t_c = xi + O(xi^2). On the other branch, log y = log(-y) - i pi, the printed formula equals t_cf - sqrt3 pi - i sqrt3 (l + g), and so does not vanish at xi = 0.

**Corollary G (refined gap, BS rmk:cfgap).** [P mod Proposition 13.2_P2, Lemma (t_cf), and draft Theorem R_P2, Corollary SV_P2, Theorem VI. The latter three are P mod (MV) = LLZ eq. (1) and textbook facts, s.1.3. Through Proposition 13.2_P2, also P mod eq:J=I (CCIT15 Thm 31, s.1.2).] For g >= 2, with GW^cf_g := Cf(F_g^{T_F}), a germ at xi = 0, and with t_cf the analytic continuation to a neighbourhood of xi = 0 of BS's expression on (-1/27,0) with log y = log(-y) + i pi, which equals t_c by Lemma (t_cf),

    GW^cf_g = (-3)^{g-1} beta_{2g-2}(eps1,eps2)/t_cf^{2g-2} + O(1).

O(1) means: the difference is holomorphic in a neighbourhood of t_cf = 0, equivalently of y = -1/27, with no negative powers of t_cf, and its Taylor coefficients lie in K. The value of the O(1) term is not determined. The statement holds with either convention for d = 0, since the constant is O(1).

**Corollary G_{k,h}.** [Same status, including eq:J=I (CCIT15 Thm 31), plus the polynomiality of Theorem FB_B(b).] For k+h >= 2,

    F^cf_{k,h} = 𝔟_{k,h} t_cf^{2-2h-2k} + O(1),  𝔟_{k,h} = (-3)^{k+h-1} [E^{2k}(eps1 eps2)^{h-1}] beta_{2k+2h-2}.

This is (Conifold asymptotics) of BS def:dirint. Here Pi_{gamma_cf} is read through its displayed expansion: BS's Pi_{gamma_cf} solves L2b o theta = 0 and vanishes linearly at xi = 0 with leading coefficient 1 (L2b o theta is eq:pfinf, whose non-constant solutions BS identifies with periods, BS 2755-2760; the expansion is BS 929; A_cf = theta Pi_{gamma_cf} is BS 937; the asymptotics are BS 955-957). Since L2b has exponents (0,0) there, the solutions vanishing at xi = 0 with no logarithm are the multiples of t_c, so Pi_{gamma_cf} = t_c = t_cf. The coefficient 𝔟_{k,h} depends on this normalisation; the integrality wording of BS 927-929 is not used. F -> F_cf in def:dirint is ev_cf on the polynomial f_j of C Step 7 (Theorem FB_B(e) after extracting monomials in eps), whose coefficients lie in Q[X^{+-1},(1-X)^{-1}] rather than BS's Q[X^{+-1}] (s.0.3). In particular F^cf_{0,h} and F^cf_{h,0} have the leading terms of thm:cfasym with remainder O(1). The limits are coefficient extractions, so E = 0 and eps1 = 0 need not be points of T_F.

**Remark (reconciliation).** [P mod PPZ Props 2.8, 2.9. Not an input.] The BS eq:betagiv graph sum equals Z[g,0] = Fpot_g(q_0). So eq:betagiv is equivalent to Theorem VI.

### 0.3 Scope and what is not proved

Scope: the fibrewise torus T_F (eps1 eps2 E != 0), per colouring. There is no extension to general CY-torus polar parts.

Not proved here:
- anything off T_F;
- the value of the regular part of GW^cf_g;
- BS's F_{k,h} in Q[X^{+-1},F]. Only F_g^{T_F} in K[y^{+-1},Delta_c^{-1}][F] is proved. Regularity at y = 0 is established for the total, and BS 4071-4076 argues from the total, not coefficient by coefficient. Stored data at genus 2-4 have BS's form [E]. Open for all g, independent of Corollary G;
- (Orbifold regularity) of BS def:dirint;
- the holomorphic anomaly equations;
- uniqueness in BS def:dirint;
- the proofs of (MV), of Teleman's Theorem (2), of Harer stability, Looijenga and Madsen-Weiss, and of eq:J=I (CCIT15 Thm 31), which stay cited;
- that the T_max use of eq:J=I (CCIT15 Thm 31) falls under ruling Q1 (item O1). This matches draft s.0.4.

Dependency: Stage A -> Stage B -> Stage C (Proposition 13.2_P2). Lemma (t_cf) is independent of A and B. Its logarithmic form gives F_cf = XX in Proposition 13.2_P2; its constant (Barnes) gives t_c = t_cf in Corollary G. Corollary G = Proposition 13.2_P2 + Lemma (t_cf) + Corollary 2.8 + Theorem R_P2 + Corollary SV_P2 + Theorem VI, with Theorem FB_B(e) for the clause on Taylor coefficients.

End-to-end inputs of Corollary G: eq:J=I (CCIT15 Thm 31, s.1.2) and the genus-0 QDE (derived in s.1.2 from Pandharipande math/9806133v1 Lemma 1 and A1.4; eq:SHb is derived from them in A2.1); Teleman Theorem (2), Prop 8.7(ii), 2.18, 2.19; PPZ Defs 2.2, 2.5, 2.13, Prop 2.14 (and Props 2.8, 2.9 only for the remark); Graber-Pandharipande; the fundamental-class and splitting axioms of the virtual class, with Coates-Givental App. 1 (iv)-(vi) for the twist (Behrend-Manin not opened); Mumford's theorem as in Faber-Pandharipande, and their Theorem 4; DLMF 15.6.7, Barnes's second lemma, Gauss's digamma values; (MV) = LLZ eq. (1) with the textbook facts of s.1.3; and, proved in the draft, SP1-SP5, Lemmas 2.1, 2.3-2.5, Theorem A, Cor 2.2, Cor 2.8, Props 2.10, 2.12, Theorem LC_P2. BFGW Lem 4.2 is not among them.

Every link has had at least one referee pass, in its report form or, for SP1-SP5 and Lemmas 2.1, 2.3-2.5, in an earlier round; s.1.4 lists exactly what those passes did not cover. The genus-0 QDE derivation of s.1.2 was written in round 6e by the report verify:ccit and read once independently by the report refQDE (round 6f; round6_scripts/round6f_result.json). The chain as assembled here had the passes recorded at the marker ROUND6C-REFEREES in s.0.4 of the draft.

---

## 1. Inputs

"Opened" means the text was read in round 6 by a prover or referee of this stage, with who stated where it matters. "Not opened" means used in standard form.

### 1.1 BS

Allowed by the rulings, and used as definitions or inputs:
- eq:Ifun (2621), eq:pf (2628), eq:Ian (2647), eq:ADEFXdef (2681), eq:etaflat (2419): definitions, read. The Picard-Fuchs recursion for I was checked by hand.
- eq:J=I (2663) and the genus-0 QDE z d S = (dt*) S, with eq:SHb (3305): read and used. BS cites eq:J=I to MR2486673 (Coates, Commun. Math. Phys. 287 (2009) 1071-1108, arXiv:0810.2200v1, opened in round 6d). That paper proves it only on T_F (Prop 4.1, Cor 4.2). At generic alpha the source of record is CCIT15 Theorem 31, s.1.2 (item O1), and the genus-0 QDE at T_max is derived in s.1.2. eq:J=I is not re-proved.
- eq:spotP2 (4301) and BS's Asym (4256; Laplace expansion 4245-4252): definitions of the oscillating integral, used only through A3.5 and B2.
- eq:Rcf, eq:F2ref (4217), rmk:cfgap (5422), lem:cfpole (5407), thm:cfasym (5301), eq:betagiv (5122): comparison only (C Step 8, C Step 9, Appendix).
- prop:algindABC (2854) and eq:XFABC (2860): used only for statements involving A (ev_A on 𝐆 with 𝖠, Cf(R_{ai}) and the legged form of eq:gtcf). They are not used in Proposition 13.2_P2 or Corollary G.

Not inputs: rmk:refmm, BKMP/FLZ remodelling, prop:cfgraphsum, thm:cfasym, lem:cfpole, TR on the mirror curve.

### 1.2 Literature

- **eq:J=I at generic alpha on T_max: source of record (item O1; source settled in round 6e; T_max use under Q1 confirmed by the user, 2026-11-02).** Coates-Corti-Iritani-Tseng (CCIT15), A mirror theorem for toric stacks, Compos. Math. 151 (2015) 1878-1912 (MR3414388), arXiv:1310.4163v5, Theorem 31 with Definition 28, eq. (19), eq. (4) and Remark 3. Read at arXiv v5 in rounds 6d and 6e. The journal print was not read, so the number Theorem 31 in print rests on a search snippet of the Imperial copy. No erratum was found on the Cambridge page or in one search; this does not exclude one. The structure of the proof (Theorem 41, Props 45-47) was read; the proof was not re-derived. *Scope.* Z = K_P2 x C^2 is the smooth toric variety with N = Z^5, rays b_0 = (1,0,1,0,0), b_1 = (0,1,1,0,0), b_2 = (-1,-1,1,0,0), b_3 = (0,0,1,0,0), b_4 = (0,0,0,1,0), b_5 = (0,0,0,0,1), and maximal cones spanned by b_3, b_4, b_5 and two of b_0, b_1, b_2. These cones are unimodular, so Box = {0}, and L = Z(1,1,1,-3,0,0). Z is semi-projective in the sense of CCIT15 s.3.1: it has three fixed points, and Z -> Spec H^0(O_Z) = C^3/mu_3 x C^2 is projective, with -D_3 = 3H relatively ample. In fan terms the support Cone(b_0,b_1,b_2) x R_{>=0}^2 is convex and full-dimensional, and the piecewise linear function with value -1 at b_3 and 0 at the other rays is the maximum of its linear pieces, strictly. NE(Z) is generated by the line. *Torus.* CCIT15's torus is N x C^* = (C^*)^6/C^*. BS's T_max (weights (alpha_0,alpha_1,alpha_2,0,-eps1,-eps2) on C^6, BS line 2349) maps onto it with kernel mu_3. The character lattices differ by index 3, so S_T = C(alpha_0,alpha_1,alpha_2,eps1,eps2) = Q_{T_max}, and the localisation formulas agree. *Dictionary.* For the T_max action the toric divisor classes are u_i = H + alpha_i (i = 0,1,2), u_3 = -3H, u_4 = -eps1, u_5 = -eps2, with H = c_1 of O(1) with its BS lift, H|P_i = -alpha_i. So H = -u_3/3, alpha_i = u_i - H, eps_j = -u_{3+j}. The non-zero u_j|P_i (j = 0,...,5) are, as a set, BS's tangent weights at P_i (BS line 2383). The indices differ: BS's w^{(i)}_i = 3 alpha_i is u_3|P_i, and BS's w^{(j)}_i = -eps_{j-2} (j = 3,4) is u_{j+1}|P_i. The relation u_0 u_1 u_2 = 0 is BS's prod(H + alpha_i) = 0. *Substitution.* Take S empty and all t_i = 0 in (19). Then lambda = d(1,1,1,-3,0,0), d >= 0, and the factors are 1/prod_{m=1}^d (H+alpha_i+mz) for i = 0,1,2, prod_{m=0}^{3d-1}(-3H-mz) for i = 3, and 1 for i = 4,5. So I_CCIT(Q,z) = z Sum_d Q^d I_d(z), where z y^{H/z} Sum_d y^d I_d(z) is eq:Ifun. *Mirror map.* Each summand with d >= 1 is -3H times a function homogeneous of degree -1 in (H,alpha,z), with expansion at z = infinity beginning 3H(-1)^d (3d-1)!/(d!)^3 z^{-1}. Reducing H^3 by the relation keeps the power of z. Hence I_CCIT(Q,-z) = -z + g(Q)H + O(1/z), with g(Q) = 3 Sum_{d>=1} (3d-1)!/(d!)^3 (-Q)^d independent of alpha and g(0) = 0; g = I_0^{[1]} - log y in the notation of eq:Ian. Theorem 31 says this is a point (4). The H_+ part of (4) is -z + t(z), so t(z) = g(Q)H, and by Remark 3 I_CCIT(Q,z) = J(g(Q)H,z). *Divisor equation.* J(tau H,z) = z e^{tau H/z}(1 + Sum_d (Q e^{tau})^d <<phi_a/(z(z-psi))>>_{0,1,d} phi^a), including the d = 0 classical terms z e^{tau H/z}. Hence Sum_d Q^d I_d = e^{g(Q)H/z}(1 + Sum_d (Q e^{g(Q)})^d <phi_a/(z(z-psi))>_{0,1,d} phi^a). Rename Q as y and multiply by z y^{H/z}. This is eq:J=I with Q = y e^{g(y)} = exp I_0^{[1]}, by eq:Jred. *Comparison with A1.2.* CCIT15 s.2.2 defines the invariants on Z by virtual localisation over the proper fixed locus (via Liu, arXiv:1107.4712v3, eq. (51)), which is BS's Notation at line 1233 for this Z and T. A1.1 and A1.2 identify this with the twisted integral on Mbar_{0,n}(P2,d) in genus 0. For d = 0, Mbar_{0,n}(Z,0) = Mbar_{0,n} x Z, and both sides are Sum_i prod a_j|P_i /W_i^0 times a psi integral on Mbar_{0,n}; the pairing is BS's eta. On genus-0 trees R pi_* f^* O_w = O_w, so the C^2 factor multiplies each invariant by 1/(eps1 eps2) and eta^{ab} by eps1 eps2, and J is unchanged. *Inputs of CCIT15's proof.* Brown's recursive characterisation of the cone in the form of Theorem 41; Liu's localisation formula; the Coates-Givental/CCIT hypergeometric modification with Tseng's quantum Riemann-Roch at the fixed points (here BG_sigma is a point, so that theory is on Mbar_{0,n}); and Theorem 2, whose TRR is the genus-0 psi-class relation. Eynard-Orantin topological recursion, remodelling and mirror curves do not occur; a text search finds no such term, and [FLT12] is cited only as an application upstream of nothing used in the proof. *Second source, read only at the level of the statements, not adopted.* Ciocan-Fontanine-Kim, arXiv:1304.7056v3, Theorem 1.2.2 and Corollaries 1.4.1, 1.4.2 (targets W//G with isolated torus-fixed points, semi-positive case). *Other sources read in round 6d and found not to cover T_max.* MR2486673 = Coates, arXiv:0810.2200v1, Prop 4.1 and Cor 4.2, prove eq:J=I on T_F only (lambda = E, p = H + E/3), since its torus covers the trivial action on P2. Coates-Givental (Ann. of Math. 165 (2007) 15-53, math/0110142) Theorem 2 and CCIT (Duke Math. J. 147 (2009) 377-438, math/0702234) Theorem 4.6 have T acting on fibres only. Givental alg-geom/9701016 Theorem 0.2 has compact X and non-negative bundles. Givental alg-geom/9603021 Theorem 9.5 gives the T-equivariant small J-function of P^n only. *Evidence [E], not a proof step.* s.A.1 records the exact check of eq:J=I through y^4 at generic alpha.
- **Genus-0 QDE at T_max [P mod the following two statements].** No numbered theorem at T_max was found: Iritani arXiv:0903.1463v3 Prop 2.4 assumes X proper, and Pandharipande math/9806133v1 s.1 assumes X convex. The derivation is short. Input (a): Pandharipande, Rational curves on hypersurfaces (after A. Givental), math/9806133v1, Lemma 1: for n >= 3, psi_1 on Mbar_{0,n}(P2,d) is the sum of the boundary divisors separating marking 1 from markings 2, 3 (P2 is convex; opened in round 6e). Lemma 1 is stated non-equivariantly; the T-equivariant form is as follows. Let c: C -> M x P^1 be the contraction of the universal curve onto the stabilisation of (C; p_1, p_2, p_3), with w the coordinate on P^1 such that p_1, p_2, p_3 go to 0, 1, infinity. T acts on M by composition with f only, so c is T-invariant and T acts trivially on P^1. Hence s_1^* c^* dw is a T-invariant section of L_1. It vanishes exactly where c contracts the component of p_1, that is on the divisors separating 1 from 2, 3. The vanishing order is 1 at the generic point of each: near a node xy = s, with p_1 at y = y_0, one has w = s/y - s/y_0 and dw = -s dy/y^2. So c_1^T(L_1) = Sum [D]^T, with no character correction. Input (b): the boundary structure with the normalisation sequence, as in A1.4 (restriction of e^{-1}(R pi_* f^* N) to a boundary divisor is the product of the two factors times e(ev^* N), and e(N)[Diagonal of P2] = Sum eta^{ef} phi_e (x) phi_f). The restriction holds for every stable splitting, because 0 -> f^*N -> f_A^*N + f_B^*N -> ev^*N -> 0 is exact on the glued family and R pi_* commutes with base change. If d_1 = 0, side A contributes e(ev^*N)^{-1}, which cancels the node factor; e(N) = -3H eps1 eps2 is invertible in V because H|P_i = -alpha_i != 0, so the TRR is an identity over Kw with no non-equivariant limit. The case d_1 = 0, A = {b} is not a divisor; the term eta(e,b) of V takes its place, through Sum eta(e,b) eta^{ef} <<f,c,a>> = <<b,c,a>>. All integrals are over the smooth proper Mbar_{0,n}(P2,d), with the fundamental class as virtual class, so no properness of Z is used. Z enters only through A1.1 and A1.2: the fixed loci, universal curves and cotangent lines coincide, so psi insertions are covered, and BS's pairing is eta. This gives the genus-0 topological recursion relation (TRR) of the theory of Z, namely <<c, a, psi^{k+1} b>> = Sum_{e,f} <<psi^k b, e>> eta^{ef} <<f, c, a>>, with only stable splittings on the right. This is the psi-class relation on Mbar_{0,n}; it has nothing to do with Eynard-Orantin TR. Here <<x_1,...,x_n>> := Sum_{d>=0} Sum_{m>=0} (Q^d/m!) <x_1,...,x_n,t,...,t>_{0,n+m,d}, the sum over (d,n+m) with d > 0 or n+m >= 3, in Kw[[Q]][[t_0,t_1,t_2]], with Q a Novikov variable. This is not eq:gwprim as printed (d > 0, no 1/m!, no Q^d). It is the reading that eq:qprod and eq:Jred use, and the summation of Pandharipande s.1.3. Put V(a,b) := eta(a,b) + <<a, b/(z-psi)>>, so that S(a) = Sum H^c eta^{cb} V(a,H^b) by eq:fundsol. Then z d_c V(a,b) = <<c,a,b>> + Sum_{k>=0} z^{-k-1} <<c,a,psi^{k+1} b>> = Sum_{e,f} V(e,b) eta^{ef} <<f,c,a>> = V(c*a, b). So z d_c S(a) = S(c*a). By the string equation (Pandharipande s.1.2 I, from psi_i = ft^*psi_i + D_{i,n+1}, twist pulled back) and <1,a,b>_{0,3,0} = eta(a,b), one has 1*a = a and S(1) = J/z, and S(H) = d_1 J by eq:fundsol. So S(H) = z theta_Q S(1) in A2.1 is the case c = H, a = 1. Restrict to t_0 = t_2 = 0. For d > 0 the descendant divisor equation <x, psi^k b, H>_d = d <x, psi^k b>_d + <x, psi^{k-1}(bH)>_d holds (from psi_2 = ft^*psi_2 + D_{2,3}, ft_* ev^*H = d, and the twist pulled back, Coates-Givental App. 1 (iv)). For d = 0 direct evaluation on Mbar_{0,2+m} gives <a, psi^{m-1} b, H,...,H>_{0,2+m,0} = eta(a H^m, b), other powers of psi vanishing. Since multiplication by H is eta-self-adjoint, S(a)|_{t_1 H} = e^{t_1 H/z}(a + Sum_{d>0} (Q e^{t_1})^d Sum_{c,b} H^c eta^{cb} <a, H^b/(z-psi)>_{0,2,d}). Each coefficient of z^{-k} lies in Kw[t_1][[Q e^{t_1}]]. So we may set the Novikov variable to 1, call Q e^{t_1} the new Q and put t_1 = log Q. Then d_1 acts as theta_Q, including on Q^{H/z}. WDVV is not used. Since z d_c [S(e*a)] = z S(d_c(e*a)) + S(c*(e*a)) and d_c d_e = d_e d_c, we get z S(d_c(e*a) - d_e(c*a)) + S(c*(e*a) - e*(c*a)) = 0. S = 1 + O(1/z) is invertible and both arguments are free of z, so c*(e*a) = e*(c*a). This uses the full formal t-space, not only the small locus. Same form of (6) and of the fundamental solution: Pandharipande s.1.2 (6), s.1.3 Prop 2. Consistent with CCIT15 Theorem 2 with (6) and Givental math/0305409v1 Theorem 1 (TRR, string and dilaton equivalent to the cone property); these were read and are not used as inputs. Evidence [E]: s.A.1.
- **Teleman, arXiv:0712.0160.** Opened (provers and referees A, A2; v3 text). Numbers confirmed: Theorem (2) of the Introduction; Definition 7.1 (shift with sign (-1)^m); opening of s.7 ("possibly formal" neighbourhood); (7.5), (7.6) (inputs E_u^{-1}(psi_i)); Proposition 7.7 with (7.7.a), (7.7.b); Remark 7.8; Example 7.17 (divisor equation; formal Frobenius manifold if convergence fails; its printed sign is not relied on, see A2.8); Proposition 8.2; Remark 8.3; Proposition 8.4; (8.6); Proposition 8.7(i),(ii); 2.17 (Harer stability), 2.18 (Looijenga), 2.19 (Madsen-Weiss). Originals pinned in round 6d. Looijenga, J. Algebraic Geom. 5 (1996) 135-150, arXiv:alg-geom/9401005v2 (opened), Props (2.1), (2.2): H^*(M_g;Q)[u_i] -> H^*(M_g^X;Q) is an isomorphism in degree <= N(g). M_g^X is the open moduli space of curves with distinct marked points, u_i = c_1 of the relative tangent sheaf = -psi_i, and N(g) >= g/2 - 1 (Ivanov, as quoted in Looijenga s.1). Madsen-Weiss, Ann. of Math. 165 (2007) 843-941, arXiv:math/0212321v3 (opened), Theorem 1.1.1 and the deduction of Mumford's conjecture after it; Teleman's item [MW] is their ECM 2005 survey, not this paper. Harer, Ann. of Math. 121 (1985) 215-249, and Ivanov, Leningrad Math. J. 1 (1990) 675-691: not opened; their ranges are as reported by Looijenga s.1 and Wahl, arXiv:1006.4476v2, s.1. Range: for integer k, k < (g-1)/2 is equivalent to k <= g/2 - 1, which is Ivanov's bound. It is conservative against Wahl Theorems 1.1, 1.2 (Randal-Williams, arXiv:0909.4278v4, s.1.3; 3k <= 2g - 2, closed case 3k <= 2g). A2.8 needs only some range growing with g, so the choice of bound plays no role. Madsen-Weiss cannot be removed from the chain: Theorem (2) of Teleman uses it (his abstract and (1.5)), although A2.8 itself uses only injectivity. **The proof of Proposition 7.7 was read, not re-derived; it is a sketch in places. Since round 6c it is not an input:** A2.8 derives the flatness equation directly, and (7.7.a) appears there only as a comparison. Theorem (2) is a non-classical citation on which the graph sum rests; eq:J=I is a further citation. Journal numbering was not compared with arXiv v3.
- **Pandharipande-Pixton-Zvonkine, J. Amer. Math. Soc. 28 (2015) 279-309, arXiv:1303.1043v3.** Opened (provers, referees A2, C, C2; numbers as in the arXiv version). Numbers confirmed: Definitions 2.2, 2.5, 2.13; Propositions 2.3, 2.8, 2.9, 2.14 (left group action of the unit-preserving action); the line on "the same convention" as Givental and Teleman; s.3.2; eq. (12); footnote 10; Remark 1.6.
- **Graber-Pandharipande, Invent. Math. 135 (1999) 487-518 (journal data recalled, not verified), arXiv:alg-geom/9708001v2.** Numbers as in the arXiv version: formula (1) in the Introduction; Appendix A (global embeddings), B (obstruction theory), C (localization for DM stacks). Provers did not open it. Referee B opened it: formula (1) and Appendices A, C are stated for a C^*-action, with the embedding hypothesis confirmed for Mbar_{g,n}(V,beta) after a cover of C^*. The passage to T is in B5.
- **Mumford.** Original opened in round 6d (OCR scan via INSPIRE): D. Mumford, Towards an enumerative geometry of the moduli space of curves, Arithmetic and Geometry II, Progr. Math. 36 (1983) 271-328, formula (5.2), Corollary (5.3) (ch_{2k} 𝔼 = 0), (5.4) with (5.5) (c(𝔼) c(𝔼^vee) = 1 on M_g, extended to Mbar_g), pp. 304-307. It treats unmarked curves, uses another Bernoulli convention, and FP s.1.1 records a typographic error in its boundary term, so the form of record stays FP's. The theorem is used in the form of Faber-Pandharipande s.1.1 (general family pi: C -> M), opened. The marked form is derived in A3.2 from that statement. Coates-Givental (21) was opened but its OCR is partly garbled, and it is no longer relied on.
- **Faber-Pandharipande, Hodge integrals and Gromov-Witten theory, math/9810173 (journal data recalled as Invent. Math. 139 (2000) 173-199, not verified).** Opened: Theorem (Mumford) in s.1.1; Theorem 4, which reads Int_{Mbar_g} lambda_{g-1}^3 = (|B_{2g}|/2g)(|B_{2g-2}|/(2g-2))/(2g-2)!, g >= 2; c(𝔼) c(𝔼^vee) = 1, quoted in s.2 just before eq. (22). With lambda_{g-1}^2 = 2 lambda_g lambda_{g-2} this gives Int lambda_{g-2} lambda_{g-1} lambda_g = half that value.
- **DLMF 15.6.7.** Prover: from memory. Opened by referees C and C2: it matches, with condition |ph(1-z)| < pi and the contour separating the poles of Gamma(a+t), Gamma(b+t) from those of Gamma(c-a-b-t), Gamma(-t).
- **Barnes's second lemma.** E. W. Barnes, Quart. J. Math. 41 (1910) 136-140: **not opened**. It is not in DLMF s.5.13 (which has only the beta integral 5.13.3). Bailey, Generalized Hypergeometric Series, Cambridge Tract 32 (1935), s.6.2, pp. 42-43: section and pages taken from MathWorld (BarnesLemma), book not opened. Andrews-Askey-Roy, Special Functions (1999), Theorem 2.4.3: number taken from arXiv:2107.13871v2 App. A (A.3), book not opened. The form in C Step 6 is the one displayed on MathWorld (beta_1 = d, beta_2 = e, beta_1 + beta_2 = alpha_1 + alpha_2 + alpha_3 + 1) and on Wikipedia (Barnes integral). The symmetric form is Jantzen, arXiv:1211.2637, eq. (2) (opened), with lambda = (1/3, 2/3), rho = (0,0,0), z = -s, 0 < Re z < 1/3. Substitution checked by hand: under s = -(sigma + 2/3) the orientation change is absorbed by ds = -dsigma, so the upward contour in s becomes the upward contour in sigma with sign +, and -2/3 < Re sigma < -1/3 corresponds to -1/3 < Re s < 0. Numerical confirmation is in s.A.3 [E]. The constant t_c = t_cf rests on this lemma and numerics are not a proof step, so before submission the form should be compared with a printed source.
- **Gauss's values of psi at 1/3, 2/3; reflection formula for psi'.** Standard, not opened.
- **BFGW, arXiv:2001.05347.** Opened by referee C: Lemma 4.2 reads "S, X, I11 are algebraically independent over ℂ", with I11 = A, S = F - (X-1)/3 (their (16), (57), (58)). Independence of A, B, C is in their s.4.1 (Lemma 4.1 plus Zagier Prop 20), and Lemma 4.2 is its corollary. They cite the identification (57) to other papers, so eq:ABCmod sits underneath it. Used only for statements involving A. Zagier and the references for (57) were not opened.
- **Coates-Givental math/0110142 App. 1; Givental math/0008067.** Opened by the prover of Stage A; used for the unit and gluing axioms (CG (iv)-(vi)) and as comparison (Givental (3)-(6), Proposition (a)-(d), s.2.3 Lemma). Givental's Theorem 1 does not apply directly, since Z is not compact. Behrend-Manin axioms: not opened.
- **Lho-Pandharipande, Stable quotients and the holomorphic anomaly equation, Adv. Math. 332 (2018) 349-402 (BS LP18:HAE), arXiv:1702.06096v4. Opened; comparison only.** They work at (19): lambda_i = zeta^i lambda_0, so e_1(alpha) = e_2(alpha) = 0, with alpha_i = -lambda_i, q = y, and the K_P2 factor only. Their (24) is eq:Ifun. Their rows (27), (30) have the form of A2.1, A2.5 at these weights: Y_i = lambda_i L, L = (1+27y)^{-1/3}; W(Y_i) = 3 e_3(alpha) is constant, so Psi_{0i} is constant, matching R_0 = 1; S(H), S(H^2) carry L lambda_i/C_1 and L^2 lambda_i^2/(C_1 C_2), with C_1 = A, C_2 = D. Their (29), C_1^2 C_2 = X, is the relation A^2 D = X of Prop 2.12. In Lemma 8 the coefficient of X_LP := D C_1/C_1 = F in R_{2k} is -R_{1,k-1}/L, which is the term -z F S(H) from S(H^2) = z theta_y S(H)/D at e_1 = 0 (A2.1), checked by hand. Their R_k vanish at q = 0 for k >= 1, with the Hodge class (36) kept at the vertex; this is R~ in R = R~ diag(h_wt) of A3.6, and (36) matches A3.1 with the C^2 weights removed. Their route is stable quotients, wall-crossing and Givental, not Teleman, and only at this subtorus; it is a cross-check of form and constants, not of the C^2 factor or of T_F. R_1 = (1-L^2)/18, R_2 = (1-24L-2L^2+25L^4)/648 were compared exactly with eq:Ifun at H = lambda_0 = 1 through y^8 [E].
- **Not used and not opened:** CCIT MR4047552 Prop 6.9, HK10, Coates-Iritani 2018, BFGW Conj 1.9. That the E = 0 and eps1 -> 0 extractions are those conjectures is BS's statement (line 5414), not checked here. The identification of F_{0,h} with the unrefined local P2 potential is in C Step 7; it is not used in the proof.

### 1.3 Draft inputs (current status, checked in the current draft s.0.2, s.7.6, s.8.1, s.8.2)

In s.1.3, s.1.4 "this round" means round 6b. Round 6a (D_rows, D_core, D_hodge, D_VI) re-derived Theorem A, Prop 2.10 S0-S6, Theorem LC_P2, s.7.5, Corollary SV_P2 and s.7.8 by hand; see draft s.0.4.

| Item | Status in the draft | Round 6 coverage relevant here |
|---|---|---|
| SP1-SP5 (s.2.2) | S, written out | Not re-proved line by line by any referee of this round. SP5 parameters, the Barnes form and sign spot-checked by hand (B, C, C2). Consequences tested through z^9 [E]. Two referees in round 3, as F1-F5 of R3 Part C. Round-4 critic 1 re-derived Lemma 2.4 and did not re-derive the F1-F4 calculus in full. No pass in the present text in round 6b. In round 6c, refS2 re-derived the calculus as used by FB:A3.5, B2, B4 (marker ROUND6C-REFEREES in s.0.4 of the draft). |
| Lemma 2.1 | P | Re-derived by hand (Stage B prover, referee B). |
| Lemmas 2.3, 2.4 | P | Lemma 2.3 re-derived by the prover of Stage B. Lemma 2.4 read and used; tested numerically [E]. |
| Lemma 2.5 | P | (i) recomputed from scratch (A); kernel argument re-read (B). |
| Theorem A, Cor 2.2 | P | Theorem A re-derived in round 5; read in this round. |
| Thm 2.7(e), steps 1-3 | P mod eq:J=I (CCIT15 Thm 31, s.1.2), genus-0 QDE (derived in s.1.2) as a whole theorem; steps 1-3 use only SP1-SP5, Lemma 2.4, Cor 2.2 and are not on the route except as FB:B4 | The stationary-phase part, with the s-chart, is now written in B4 and re-derived by hand (B, C, C2). (C1)-(C3) are replaced by Stage A. |
| Prop 2.10 | P | Referee passes in rounds 5 and 6a. S0 and (Pair) re-derived (B, C, C2). S1-S6 not re-derived in this round. |
| Prop 2.12 | P, no citation | Round 5 and 6 referee passes. |
| Cor 2.8 | P | Algebra on the closed form. |
| Theorem LC_P2 | P | Round 6 pass (D_core). |
| Lemma A_P2(b) | **P mod (MV)** | lemmaA/LEMMA_A_PROOF.md Theorem L; referee L_logic, sound with minor fixes. |
| Theorem R_P2 | **P mod (MV)** | Through Lemma A_P2(b) for 2 <= h <= g-1. Not re-derived in this round; per-colouring bounds tested at g = 2, 3 [E]. |
| Corollary SV_P2 | **P mod (MV)** | Through Lemma A_P2(b) at h = g and Theorem VI. Tested, not re-derived. |
| Theorem VI | **P mod (MV)** | LLZ eq. (1) only, at mu = 1^a; LEMMA_A_PROOF.md Prop 2.1, Lemma 3.2, (4.2) at k = 2; textbook facts (string, dilaton, Frobenius, Murnaghan-Nakayama, principal specialisation); (7.2) of the draft. Two referee passes (rounds 5, 6). |

**(MV)** is the Marino-Vafa formula, LLZ eq. (1) (eqs. (1)-(3) of LLZ, J. Differential Geom. 65 (2003) 289-340, p. 291, numbered as in arXiv:math/0306434v3; JDG print opened in round 6d through the EMIS mirror: (1)-(3) on printed p. 291, and the name Marino-Vafa formula for (1) on p. 292. LLZ Theorems 1, 2 are the cut-and-join equations, not (MV). Okounkov-Pandharipande: Geom. Topol. 8 (2004) 675-699, Theorem 1, opened), with Okounkov-Pandharipande Theorem 1 as a second published proof. Its proof was not re-checked. Textbook facts for Lemma A: string and dilaton equations for pulled-back classes, the genus-0 multinomial formula, Frobenius formula, Murnaghan-Nakayama, the hook-length form of the principal specialisation, Lagrange inversion. **LLZ eq. (5) and the boson-fermion correspondence are no longer inputs** (the report FB_C still lists them; that list is out of date). The first mKP equation (4.2) is derived in LEMMA_A_PROOF.md s.4 from Murnaghan-Nakayama and Frobenius.

### 1.4 What the round-6b passes did not cover

- Proofs of Teleman's Theorem (2) and of 2.17-2.19 (Harer stability, Looijenga, Madsen-Weiss): used as stated, not re-derived. Prop 7.7 is not an input since round 6c; its proof was read, and its use is now a remark in A2.8.
- SP1-SP4 line by line, and the proof of Lemma 2.4 (see the table).
- eq:J=I itself; GP's proof; (MV).
- Prop 2.10 steps S1-S6 and Theorem A in this round; Theorems R_P2, SV_P2, VI (tested, not re-derived, in this round).
- Referee B's verdict is conditional on Stage A, and referee C2's verdict on Stages A and B. Stage A's former formal-family use of Prop 7.7 was examined by referees A, A2 only. The direct derivation now in A2.8 has had two independent hand re-derivations, recorded in A2.8. Referee C also took Theorems LC_P2, R_P2, SV_P2, VI as proved and read Stages A, B only at their interface.
- The class-level identity Omega = R.omega rests on Teleman. Evidence only: the round-6d computation of s.A.1 pairs it with psi classes pulled back from Mbar_{g,n} and with powers of lambda_1 at d <= 3 (g <= 2) and d <= 2 (g = 3, n = 1, psi only). It was not paired with boundary strata or with general tautological classes, and it is not a proof step.
- The assembled text. (Later passes: marker ROUND6C-REFEREES in s.0.4 of the draft.)

---

## 2. Stage A: proof of Theorem FB_A

No step uses a finite-order computation.

### A1. The CohFT

**A1.1 Fixed loci [P].** A compact irreducible T-invariant curve in Z maps to a T-fixed point of C^2, that is, the origin. A compact curve in K_P2 lies in the zero section, because K_P2 -> C^3/mu_3 is affine away from P2. So Mbar_{g,n}(Z,d)^T = Mbar_{g,n}(P2,d)^T, which is closed in the proper Mbar(P2,d).

**A1.2 All N-weights move, g >= 0, rank 3 [P].** Fix a fixed map with vertices C_v at P_{i(v)} and edges C_e covering a line P_iP_j with degree d_e. By the normalisation sequence, R pi_* f^* N|_F is:
- for O(-3): vertices (O - 𝔼_v^vee) 3alpha_i, flags -3alpha_i, and edges H^1(P^1, O(-3 d_e)) (H^0 = 0) with weights 3alpha_i + l (alpha_j - alpha_i)/d_e = ((3d_e - l) alpha_i + l alpha_j)/d_e, l = 1,...,3d_e - 1;
- for O_{-eps}: (O - 𝔼^vee)(-eps).

Every weight is a non-zero linear form. They all remain non-zero on T_F, where they become E, -eps1, -eps2. For d > 0, H^0(C, f^*O(-3)) = 0 by connectedness, and R^1 has rank 3d + g - 1. Hence the fixed part of the Z-obstruction theory on F is that of P2, the virtual classes [F]^vir coincide, and e(N^vir_Z) = e(N^vir_P2) e_T(R pi_* f^* N). By Graber-Pandharipande on the proper Mbar(P2,d) (with the passage from C^* to T in B5), BS's localisation definition (Notation at line 1233) equals the twisted integral in the statement.

**A1.3 Regularity on the diagonal [P mod GP; this is B5].** The report's series argument through Q[mu^{+-1},eps^{+-1}][[beta]] is replaced by the polynomiality of B5, whose proof uses only A1.2 and Graber-Pandharipande, so there is no circularity. For d >= 1, eps1 eps2 GW^T_{g,d} lies in Q[alpha,eps]. For d = 0 and g >= 2, 3 e_3 eps1 eps2 GW^T_{g,0} lies in Q[alpha,eps]. Since eps1, eps2 and e_3 are not in 𝔮, GW^T_{g,d} lies in O_w, and its value at alpha* is the T_F integral. This proves BS line 4056 and (Contravariance) for this case.

(The report's route has two defects, recorded for traceability. Q[mu^{+-1},eps^{+-1}][[beta]] is not an algebra over the local ring Q[beta,mu,eps]_{(beta)}, since 1/(mu+eps) is not in it; the correct setting would be the completion Q(mu,eps)[[beta]] with bO^ ∩ O = bO. And that the rational function and the beta-series are the same element needs an argument. Neither is needed now.)

**A1.4 Axioms [P mod C].** Over Kw[[Q]]:
- S_n symmetry is clear.
- Unit: Omega_{g,n+1}(gamma, 1) = p^* Omega_{g,n}(gamma), from the fundamental-class axiom of the virtual class with E_{g,n+1,d} = ft^* E_{g,n,d} (Coates-Givental App. 1 (iv)). Omega_{0,3}(a,b,1) = eta(a,b), since d > 0 pulls back from Mbar_{0,2}(P2,d).
- Gluing: the splitting axiom plus iota^* R pi_* f^* N = (two sides, or one side) - ev^* N (Coates-Givental App. 1 (v), (vi); directly, the normalisation sequence). The node factor e(N) turns the diagonal of P2 into the bivector eta^{-1}, which is Sum_i W_i^0 e_i(0) (x) e_i(0).

**A1.5 Divisor [P mod the fundamental-class axiom].** ft_*(ev_{n+1}^* H) = d (equivariant degree-0 class, the fibre degree of f^*O(1)), so p_* Omega_{g,n+1,d}(gamma, H) = d Omega_{g,n,d}(gamma), also at d = 0. Iterating, Sum_m (t^m/m!) p^m_* Omega_{g,n+m}(gamma, H^m) = sigma(Omega_{g,n}(gamma)), with sigma := e^{t theta_Q}: K~w((Q)) -> K~w((Q))[[t]] the ring map Q -> Q e^t.

**A1.6 Semisimplicity [P].** The matrix of H* reduces mod Q to H cup, with distinct eigenvalues -alpha_i. By Hensel, the eigenvalues lambda_i and idempotents e_i lie in Kw[[Q]], and Delta_i(0) = W_i^0.

**A1.7 Freeness of the action [P].** Suppose R.omega = omega over a field of characteristic 0, with R = 1 + R_m z^m + ... in the f_i basis.
- Evaluate (R.omega)_{0,m+3}(f_i, f_j, ..., f_j), i != j. Graphs with an edge separate f_i, and translations vanish on mixed insertions, so only the leg term at the f_i-leg survives. This gives -(R_m)_{ji} Delta_j^{(m+1)/2} psi_1^m, and psi_1^m != 0 on Mbar_{0,m+3}. Hence R_m is diagonal, R_m = diag(r_i).
- Symplecticity makes m odd, m = 2k - 1 with k >= 1: R(z) R^T(-z) = 1 gives R_m + (-1)^m R_m^T = 0, so a diagonal R_m with m even vanishes.
- **Colour isolation.** Work on Mbar_{k,1} with the insertion e_i (on Mbar_k there is no marking, and Omega_{k,0} is a sum over colours, so one colour cannot be separated). The degree-m part of (R.omega)_{k,1}(e_i) is linear in R_m, monochromatic (a half-edge at a vertex carrying e_i or joined to it must have colour i), and equals Delta_i^{k-1} times the degree-m part of (exp(r_i z^m).1)_{k,1}. By A3.2 this is Delta_i^{k-1} r_i ((m+1)!/B_{m+1}) ch_m(𝔼). By the flat unit this class on Mbar_{k,1} is p^* of the class on Mbar_k for k >= 2. p^*: H^*(Mbar_k;Q) -> H^*(Mbar_{k,1};Q) is injective because p_*(psi_1 p^* x) = (2k-2) x. So ch_{2k-1}(𝔼) != 0 on Mbar_{k,1} by the first paragraph of A3.3, which uses only c(𝔼)c(𝔼^vee) = 1 and FP Theorem 4, not A1.7, A1.8 or A3.1 (for k = 1, Mbar_{1,1} is used directly). Hence r_i = 0.

**A1.8 Teleman [C, P].** F_0 := K~w((Q)) has characteristic 0 and cardinality of the continuum, so it embeds in ℂ. Since H^*(Mbar_{g,n}) is finite dimensional, iota(Omega) is a CohFT with flat unit in H^*(Mbar_{g,n};Q) (x) ℂ, semisimple because Omega is semisimple over F_0. PPZ Defs 2.2, 2.5, 2.13 and Prop 2.14 are identities of graph sums, and their proofs hold over any field of characteristic 0. Teleman Theorem (2) and Prop 8.7(ii), in PPZ's form (s.3.2), give Omega = R.omega with R over ℂ. Automorphisms of ℂ over F_0 fix Omega and omega, so by A1.7 they fix R. The fixed field of Aut(ℂ/F_0) is F_0, for F_0 of characteristic 0 (an algebraic x not in F_0 has a different conjugate, a transcendental x can be moved to x + 1, and these extend to ℂ). So R has coefficients in K~w((Q)).

### A2. The R-matrix

**A2.1 Rows are formal [P mod eq:J=I (CCIT15 Thm 31, s.1.2) and the genus-0 QDE (derived in s.1.2)].** From eq:J=I, S(1) = I/z, S(H) = theta_y I/A, S(H^2) = (z theta_y S(H) - e_1 E(y) S(H))/D (eq:SHb), and z theta_Q S(H^a) = Sum C^c_{1a} S(H^c). Further

    I|_{P_j} = z y^{-alpha_j/z} Sum_d y^d prod_{m=0}^{3d-1}(3alpha_j - mz)/prod_l prod_{m=1}^d (alpha_l - alpha_j + mz),

so s_{aj} = y^{-alpha_j/z} s^_{aj} with s^_{aj} in K~w(z)[[y]]. The d-th Q-coefficients of S(1), S(H), S(H^2) are built from I by theta and by division by A and D, so they are rational in z. eq:J=I is an identity of 1/z-series, and a 1/z-series identity between rational functions is a rational identity. This covers the third QDE row used in A2.4. So all these identities may be expanded at z = 0.

*Row 2, with no lem:modrel.* With s = 1/z, eq:Ifun gives S(1) = I/z = 1 + s H I_0^{[1]} + s^2 (H I_1^{[1]} + H^2 I_1^{[2]}) + O(s^3) (draft Prop 2.12, Step 0). So S(H) = z theta_Q S(1) = H + s (e_1 (G/A) H + (P_I/A) H^2) + O(s^2), with G, P_I as in that Step 0. Since S(H^c) = H^c + O(1/z), comparing the z^0 terms of z theta_Q S(H) = C^0_11 S(1) + C^1_11 S(H) + C^2_11 S(H^2) gives C^0_11 = 0, C^1_11 = e_1 E(y)/A, C^2_11 = D/A, by eq:ADEFXdef alone. D(0) = 1, so S(H^2) = (z theta_y S(H) - e_1 E(y) S(H))/D. Draft Remark 2.13(a) records the same values as an audit of BS; this paragraph is the proof and does not use the remark.

**A2.2 [P].** Here f_j^* must be the constant covector at Q = 0. BS line 3537 says "at t"; with a moving covector (Flatness) fails.

**A2.3 Formal WKB [P, rank-free].** Let Psi in GL_3(K~w[[Q]]) have columns the eta^{-1}-orthonormal eigen-covectors, with Psi_{aj}(0) = (-alpha_j)^a (W_j^0)^{-1/2}. W_Psi := Psi^{-1} theta_Q Psi is antisymmetric with zero diagonal, and lies in Q Mat_3 because theta_Q kills constants. Solve (z theta_Q + lambda_j) v = C v with v = Psi x, x = Sum z^n x^{(n)}:
- (lambda_l - lambda_j) x_l^{(n)} = theta_Q x_l^{(n-1)} + (W_Psi x^{(n-1)})_l;
- for l != j the factor is a unit of K~w[[Q]], so x_l^{(n)} is determined;
- for l = j: theta_Q x_j^{(n-1)} = -(W_Psi x^{(n-1)})_j lies in Q K~w[[Q]], which is in the image of theta_Q. So a solution exists with x^{(n)}(Q=0) = delta_{n0} e_j and x^{(0)} = e_j.

**A2.4 Comparison with the S-rows [P].** First define u_j. Since lambda_j(0) = -alpha_j and A(0) = 1, lambda_j + alpha_j/A lies in Q K~w[[Q]], and U_j in Q K~w[[Q]] is defined by theta_Q U_j = lambda_j + alpha_j/A. Put u_j := U_j - alpha_j log y, so that theta_Q u_j = lambda_j. (A2.5 identifies theta_y u_j = Y_j.)

Put w^ := exp(-(u_j + alpha_j log y)/z) s^_{.j} = exp(-U_j/z) s^_{.j}, in K~w((z))[[Q]]. Here log y, not log Q, is the correct logarithm, because s_{aj} = y^{-alpha_j/z} s^_{aj}. Then w^ = e^{-u_j/z} s_{.j}, and from the QDE z theta_Q w^ = (C - lambda_j) w^, with w^|_{Q=0} = Psi_{.j}(0). In the Q-expansion, the operator on the d-th coefficient has diagonal entries zd + alpha_l - alpha_j, all invertible in K~w((z)) for d >= 1, including l = j, where zd != 0 in K~w((z)). So w^ is unique and w^ = v. This proves R~ in K~w[[y]][[z]], R~|_{z=0} = Psi, and R~|_{y=0} = Psi|_{y=0}. No oscillating integral and no analytic asymptotics are used, so (C1)'s eq:asymp is replaced by this argument.

**A2.5 Y_j and Delta_j derived, with no lem:modrel [P].** Conjugate eq:pf (satisfied by I; direct check on the recursion) by e^{u_j/z}, with theta_y u_j = A lambda_j:
- at order z^0, Lhat(A lambda_j) = 0, so Y_j := A lambda_j is the root with Y_j(0) = -alpha_j, and theta_y u_j = Y_j;
- at order z^1, (theta + (1/2) theta log W(Y_j)) Psi_{0j} = 0 (draft Lemma 2.5(i)). Since Psi_{0j} = Delta_j^{-1/2} and Delta_j(0) = W_j^0 = 3 eps1 eps2 W(-alpha_j), we get Delta_j = 3 eps1 eps2 W(Y_j).

**A2.6 Unitarity [P].** G_ab := Sum_k s_{ak}(z) s_{bk}(-z) lies in K~w(z)[[Q]] (the factors y^{+-alpha_k/z} cancel), and z theta_Q G = CG - GC^T. eta also solves this (self-adjointness), and G(0) = eta. The difference vanishes by the same invertibility as A2.4. Write s = s^ diag(y^{-alpha_k/z}). Then G = s^(z) s^(-z)^T with s^ in GL_3(K~w(z)[[Q]]), since s^(Q=0) is invertible. Inverting G = eta gives s^(-z)^T eta^{-1} s^(z) = 1. Conjugating by exp(U_k/z), that is R~_{ai}(-z) = e^{U_i/z} s^_{ai}(-z) and R~_{bj}(z) = e^{-U_j/z} s^_{bj}(z), yields Sum R~_{ai}(-z) eta^{ab} R~_{bj}(z) = delta_ij. The row formulas follow by conjugating z theta_Q S(H^a) with e^{u/z}, using C^0_{11} = eta_{11}/eta_{00} = 0.

**A2.7 Dictionary [P, with PPZ read].** eta(H^a, R(z) f_j) = R_{aj}(z) is equivalent to f_j^* R^{-1}(z) H^a = R_{aj}(-z), by symplecticity. Under it, PPZ Defs 2.2 and 2.13 become exactly BS eq:leafterms, eq:edgeterm, eq:gtH.

**A2.8 Flatness of R [P mod Teleman Theorem (2) through A1.8, A1.5 (fundamental-class axiom), PPZ Defs 2.2, 2.13 and the symplectic group (5), Mumford and FP Thm 4 through A1.7 in A1.8, and Teleman 2.17-2.19 (Harer stability with Looijenga; Madsen-Weiss; originals as in s.1.2)].**

Rewritten in round 6c. This derivation does not use Teleman's Prop 7.7. It has had two independent hand re-derivations (the proposer's and that of the round-6d referee ref28hand, with PPZ Defs 2.2, 2.13 and Teleman 2.15-2.19, 7.1, 7.6-7.8, 8.2 compared against the arXiv text). The round-6d reads refNew, coldF and coldB also re-derived it by hand and found only the notation corrected here. The record of these AI passes is the marker ROUND6C-REFEREES in s.0.4 of the draft. The finite-order tests of s.A.1 are evidence, not passes. In A2.8 only, D denotes the diagonal matrix below, not BS's D(y).

*Set-up.* Work over F_0 with the derivation theta := theta_Q, using Omega = R.omega (A1.8) and A1.5. Put B(psi) := R^{-1}(psi) and a^i(psi;x) := f_i^*(B(psi)x) for x in V over Kw, so theta x = 0. On M_{g,n} (smooth curves) only the edgeless graph survives:

    Omega_{g,n}(x_1,...,x_n) = Sum_i D_i N_i(kappa) prod_l a^i(psi_l;x_l),

with D_i := Delta_i^{(2g-2+n)/2} and N_i(kappa) := Sum_m (1/m!) p_{m*} prod_k (Delta_i^{1/2} f_i^* T(psi_{n+k})), a power series in the kappa_j, j >= 1 (finite in each degree, and given by exp(Sum_j s^i_j kappa_j) with exp(-Sum_j s^i_j psi^j) = b_i(psi) := Delta_i^{1/2} a^i(psi;1), by T(z) = z(1 - R^{-1}(z)1)), with constant term 1, the same series for every (g,n).

*Pushforward.* Let p: C -> M_{g,n} be the restriction of Mbar_{g,n+1} -> Mbar_{g,n}, with sections delta_l. On C: psi_l = p^*psi_l + delta_l, psi_{n+1} delta_l = 0, psi_l delta_l = 0 (hence psi_l^k = (p^*psi_l)^k + delta_l (p^*psi_l)^{k-1}, as in Teleman (2.15)), delta_l delta_l' = 0 for l != l', kappa_j = p^*kappa_j + psi_{n+1}^j, and the only graphs with an edge are the bubbles (0,3) carrying legs l, n+1 (a point of C has at most one node, so graphs with two or more edges miss C, and the other one-edge graphs lie over the boundary of Mbar_{g,n}). The bubble vertex is omega_{0,3} with no translation leg, since p_{m*} of a class of degree >= 2m on Mbar_{0,3+m} vanishes on Mbar_{0,3} for m >= 1. Its leg factors are R^{-1}(0) = 1, since psi vanishes on Mbar_{0,3}. On delta_l the classes kappa_j restrict to p^*kappa_j and a^i(psi_{n+1};H) restricts to f_i^*H, since psi_{n+1} delta_l = 0. So B(psi_l) = B(p^*psi_l) + delta_l (B(p^*psi_l) - 1)/p^*psi_l. After p_*, using p_*(delta_l p^*c) = c, the factor becomes (B(psi_l) - 1)/psi_l with psi_l on M_{g,n}. The bubble puts (1 - B(psi_l))/psi_l applied to H*x_l on leg l of the genus-g vertex (its edge factor at psi_h = 0 is (eta^{-1} - (1 (x) B(psi_{h'})) eta^{-1})/psi_{h'}). Since Delta_i^{(2g-1+n)/2} a^i b^i = Delta_i^{(2g-2+n)/2} f_i^*(a*b), the two delta_l terms add up to (1/psi_l)[H*,B(psi_l)] x_l on leg l. Hence

    p_* Omega_{g,n+1}(x,H) = Sum_i D_i N_i Sum_l f_i^*((1/psi_l)[H*,B(psi_l)]x_l) prod_{l' != l} a^i(psi_l';x_l') + Sum_i G_i(kappa) prod_l a^i(psi_l;x_l),

with G_i = D_i N_i(kappa) p_*[Delta_i^{1/2} a^i(psi_{n+1};H)/b_i(psi_{n+1})], a series in the kappa_j, j >= 1, with coefficients in F_0 and with the number kappa_0 = 2g-2+n among them. G_i contains no psi_l. Here N_i(p^*kappa + psi^j) = N_i(p^*kappa)/b_i(psi) is used.

*The defect M.* Define M(psi) in Mat_3(F_0[[psi]]) by

    theta a^i(psi;x) - f_i^*((1/psi)[H*,B(psi)]x) = Sum_k M(psi)^i_k a^k(psi;x),

which is possible because x -> a(psi;x) is invertible. A1.5 at n = 2, where D_i = Delta_i^g, reads

    Sum_i D_i N_i [(M(psi_1)a_1)^i a_2^i + a_1^i (M(psi_2)a_2)^i] = Sum_i Gtilde_i(kappa) a_1^i a_2^i,   Gtilde_i := G_i - theta(D_i N_i).

*Stability.* By Teleman 2.17-2.19, in cohomological degree < (g-1)/2 (deg psi_l = 2, deg kappa_j = 2j), H^*(M_{g,2};Q), the rational cohomology of the open moduli space, is the polynomial ring on psi_1, psi_2 and the pull-backs of the kappa_j of M_g. On M_{g,2}, where delta_{12} = 0, the classes used here, which are the Arbarello-Cornalba classes, satisfy kappa_j = pi^*kappa_j + psi_1^j + psi_2^j (Teleman (2.15)). This is a triangular degree-preserving change of generators, so monomials in kappa_j (j >= 1), psi_1, psi_2 form a basis in that range. Over F_0 they form a basis because H^*(M_{g,2};F_0) = H^*(M_{g,2};Q) (x)_Q F_0. theta acts as 1 (x) theta, killing the monomials. All coefficient extractions below are in this basis, and kappa^0 means no kappa_j with j >= 1. Only independence in some range growing with g is used. Take the coefficient of kappa^0 psi_2^0, where N_i = 1 and a_2 = (f_i^* x_2)_i is arbitrary, and let x_1 be arbitrary. With D = diag(D_i):

    D M(psi) + M_0^T D = Gtilde(0), diagonal, modulo that degree.

Since g is arbitrary and M does not depend on g, M_n = 0 for n >= 1 (take g > 4n + 1, so that psi_1^n has degree 2n < (g-1)/2). Further (M_0)_ik = -(W_Psi)_ik + ([Lambda, R_1])_ik in the f-frame, and both terms are skew (R_1 is symmetric by symplecticity). So (Delta_i^g - Delta_k^g)(M_0)_ik = 0, and Delta_i^g != Delta_k^g because their values at Q = 0 are (W_i^0)^g, (W_k^0)^g, and, with {i,k,m} = {0,1,2}, W_i^0/W_k^0 = -alpha_i (alpha_m - alpha_i)/(alpha_k (alpha_m - alpha_k)) is a non-constant element of Q(alpha), hence transcendental over Q and not a root of unity. Any g >= 1 serves here, since only degree 0 is used. Hence M = 0, that is

    theta(Psi^{-1} B(psi)) = Psi^{-1} (1/psi)[H*,B(psi)].

*The equation.* Taking eta-adjoints, with B(psi)^dagger = R(-psi) (this is R^{-1}(z)^dagger = R(-z), the symplectic condition of PPZ (5)), Psi^dagger = Psi^{-1}, (H*)^dagger = H*, and theta eta = 0, so that theta commutes with dagger, and putting z = -psi gives

    z theta_Q(R(z) Psi) = (H* R(z) - R(z) H*) Psi.

This is an equation for R with legs R^{-1}(psi), not for R^*. With H* Psi = Psi Lambda and theta_Q U = Lambda, it is equivalent to z theta_Q S = H* S for S = R Psi e^{U/z}, which is Givental (3).

*Remark: the kappa terms carry no further constraint [P; a consistency statement, not a repair; re-derived by hand in round 6d].* With M = 0 and c_i(psi) := Delta_i^{1/2} a^i(psi;H)/b_i(psi) = lambda_i + c_{i,1} psi + ..., one has c_{i,1} = Delta_i^{1/2} Sum_k (W_Psi)_ik Delta_k^{-1/2} = theta log Delta_i^{1/2}, which matches the kappa_0 term of theta log D_i. The flatness equation applied to x = 1 gives -psi theta b_i = c_i b_i - (lambda_i + c_{i,1} psi) b_i, which is theta s^i_j = c_{i,j+1}. So Gtilde_i(kappa) = 0 identically: M = 0 implies A1.5 on M_{g,n} in full. Compare (7.7.c) of Teleman, which was not re-derived in his notation.

*Remark: comparison with Teleman's Prop 7.7 [C; not a proof step].* Put E_t := sigma(R), over F_0[[t]]. Applying sigma to coefficients gives sigma(R).sigma(omega) = sigma(Omega). Apply the embedding iota: F_0 -> ℂ of A1.8 coefficientwise in t. Then iota(sigma(Omega)) is Def 7.1 for Z = iota(Omega) and u = -tH over ℂ[[t]], iota(sigma(R)) is its data, and the coefficient of t^0 in (7.7.a) is iota of the identity over F_0 stated above. theta_Q is never extended to ℂ. Only the direction H is used, because the integrand p_* Omega_{g,n+1}(x,H) belongs to the same family. Def 7.1 has the factor (-1)^m, and the divisor equation used is the true one, A1.5; the sign printed in Teleman Ex. 7.17 is not relied on. By A1.7, applied over the field ℂ((t)), E_t is the unique such data. (7.7.a) along the line reads d_t(𝐄 Psi) Psi^{-1} = -[𝐄, H*_t]/z, and at t = 0 this is the equation above, which is also (7.7.a) along v = -H. The legs of (7.6) carry 𝐄_u^{-1}(psi_i), and the legs of PPZ Def 2.2 carry R^{-1}(psi_i), so in this comparison Teleman's 𝐄 is PPZ's R, with the same z; this is a reading of (7.6) against Def 2.2. Teleman's Prop 8.2(ii) is consistent with this. Theorem (2) is a non-classical citation from Teleman on which the graph sum rests; eq:J=I is a further citation.

*Sign check, first order only [P].* The degree-1 part of (R.omega)_{0,4} from PPZ Def 2.2 has coefficient -Delta_i r_ij on a_j b_i c_i d_i and Delta_i^{1/2} Delta_j^{1/2} r_ij on a_j d_j b_i c_i. This is exactly the derivative of C_abc = Sum Delta_i^{1/2} a_i b_i c_i with W_Psi = [dU, R_1]. This separates R(z) from R(-z). **It does not separate R from R^*(z) := R^{-1}(-z)**, because R_1 is symmetric and R^* has the same R_1 (R^*_2 = R_1^2 - R_2). So it is not an independent check of the full flatness statement.

*Remark: convention check through PPZ eq. (12) [P, hand derivation; not a proof step, since the line in Q has no Euler field].* On a conformal semisimple Frobenius manifold with Euler field E_ul, flat basis d_a of charge q_a, mu d_a = (q_a - d/2) d_a, xi = E_ul*, let Psi be the matrix of the normalised idempotents f_i in the flat basis, and bold-R := Psi^{-1} R Psi. Then S = Psi (bold-R) e^{U/z} and z dS = (dt*) S give z d_v (bold-R) + z (Psi^{-1} d_v Psi)(bold-R) = [d_v U, bold-R]. Take v = E_ul: d_v U = U; d_v bold-R_m = -m bold-R_m; and d_v Psi = -mu Psi, which follows from [E_ul, f_i] = (d/2 - 1) f_i and [E_ul, d_a] = (q_a - 1) d_a. At z^{m+1}, and after conjugating by Psi,

    [R_{m+1}, xi] = (m + mu) R_m,

which is PPZ eq. (12) at all m. For R^* the same steps give [R^*_{m+1}, xi] = -m R^*_m + R^*_m mu, a different equation from m = 1 on. So PPZ's own R is the one with S = Psi R e^{U/z}, not R^*.

*Attribution.* PPZ footnote 10 and Remark 1.6 attribute the opposite commutator sign to Teleman's opposite Euler field, not to the (-1)^m of his Def 7.1. Both amount to u -> -u, and they are two separate statements. In the proof the sign comes directly from PPZ Def 2.2, A1.5 and the adjoint step. Def 7.1 enters only in the comparison remark, where it carries the sign of that route.

**A2.9 Ambiguity [P, no citation].** X_0 := Psi^{-1} R~^{-1} R Psi satisfies z theta_Q X_0 = [Lambda, X_0]. Here R~ is the operator with matrix R~_{aj} = eta(H^a, R~ f_j), that is, X_0 is the matrix R~_{aj}^{-1} R_{aj} with R_{aj} = eta(H^a, R(z) f_j) (A2.7). Both matrices solve z theta_Q Y = C Y - Y Lambda, by A2.4 and A2.8. With X_0 = 1 + Sum X_n z^n: (lambda_i - lambda_j)(X_n)_ij = (theta_Q X_{n-1})_ij, and lambda_i - lambda_j is invertible. By induction X_n is diagonal and theta_Q-constant; the constants of theta_Q on K~w((Q)) are K~w. Symplecticity gives X_0 = diag exp(odd series). So R_{aj} = R~_{aj} exp Sum_k a_k^{(j)} z^{2k-1}, with a in K~w, and R has coefficients in K~w[[Q]].

### A3. The constants

**A3.1 Reduction mod Q [P, GP].** R (A2.9), R^{-1} (symplecticity), e_i, Delta_i^{+-1/2} and eta^{-1} are regular at Q = 0 (Delta_i(0) = W_i^0 != 0), and the graph sum is polynomial in them, so Omega|_{Q=0} = R|_{Q=0}.omega|_{Q=0}. By A2.4 (R~|_{Q=0} = Psi|_{Q=0}) and A2.9, R|_{Q=0} = diag(s_i) on the f_i(0), s_i := exp Sum a_k^{(i)} z^{2k-1}. Edges are then diagonal and graphs are monochromatic. Using 1 = Sum e_i and Sum_v (g_v - 1) + |E| = g - 1, for n >= 1 the (e_i,...,e_i) component is (W_i^0)^{g-1} (s_i.1)_{g,n}, and mixed insertions vanish. On the geometric side, d = 0: fixed loci are Mbar_{g,n} x {P_i}, all of T_{P_i}Z moves, [F]^vir = [Mbar_{g,n}], N^vir = T_{P_i}Z - 𝔼^vee (x) T_{P_i}Z, and e(𝔼^vee (x) C_w) = Lambda^vee(w). So Omega(e_i^n)|_{Q=0} = (W_i^0)^{-1} prod_w Lambda^vee(w) = (W_i^0)^{g-1} Hodge_{wt_i}. Hence s_i.1 = Hodge_{wt_i} for all stable (g,n) with n >= 1.

**A3.2 Hodge_v = h_v.1, every exponent signed [P mod Mumford].**
- *Mumford with markings.* Faber-Pandharipande's theorem applies to the family pi: Mbar_{g,n+1} -> Mbar_{g,n}. With K := c1(omega_pi), D = Sum D_i, K.D_i = psi_i D_i, D_i^2 = -psi_i D_i and D_i D_j = 0, one has (K+D_i)^{2l} = K^{2l} + psi_i^{2l-1} D_i Sum_{r=1}^{2l} (-1)^{r-1} binom(2l,r) = K^{2l} + psi_i^{2l-1} D_i, so pi_*(c1(omega(D))^{2l}) = pi_*(c1(omega)^{2l}) + Sum_i psi_i^{2l-1}. Hence, in Arbarello-Cornalba conventions kappa_{2l-1} = p_*(psi_{n+1}^{2l}),

      ch_{2l-1}(𝔼) = (B_{2l}/(2l)!) [kappa_{2l-1} - Sum_i psi_i^{2l-1} + (1/2) xi_*(Sum_{a+a'=2l-2} (-1)^{a'} psi^a psi'^{a'})],

  for all stable (g,n), including g = 0, where both sides vanish.
- c(𝔼) c(𝔼^vee) = 1 gives ch_{2k}(𝔼) = 0 for k >= 1, and Lambda^vee(v)/v^g = c_{-1/v}(𝔼) = exp(-Sum_k (2k-2)! ch_{2k-1}(𝔼) v^{1-2k}).
- Claim: for scalar h = Sum b_k z^{2k-1}, exp(h).1 = exp(X_h) with X_h := Sum b_k ((2k)!/B_{2k}) ch_{2k-1}(𝔼).
- Proof: by PPZ Prop 2.14, d/ds (e^{sh}.1) = D_h(e^{sh}.1), where from Defs 2.2, 2.13

      D_h Omega_{g,n} = Sum_k b_k [-Sum_i psi_i^{2k-1} Omega_{g,n} + p_*(psi_{n+1}^{2k} Omega_{g,n+1}) + (1/2) xi_*(Sum_{a+a'=2k-2} (-1)^{a'} psi^a psi'^{a'} Omega')].

  The legs give -psi^{2k-1} (from R^{-1}), T(z) = b_k z^{2k} gives the kappa term, and the edge gives (z^{2k-1} + w^{2k-1})/(z + w). On exp(sX_h), which is pulled back under p and additive under xi^* for k >= 1 (loops included), the projection formula gives D_h exp(sX_h) = exp(sX_h) Sum_k b_k [kappa_{2k-1} - Sum psi_i^{2k-1} + (1/2) xi_*(...)] = X_h exp(sX_h), by Mumford's formula. So exp(sX_h) solves the same graded ODE with the same initial value.
- Comparing: h_v has b_k = -(B_{2k}/(2k(2k-1))) p_{2k-1}(1/v), so X_{h_v} = -Sum (2k-2)! ch_{2k-1}(𝔼) p_{2k-1}(1/v) = log Hodge_v. This agrees with Teleman Prop 8.7(i) and Givental's Lemma in s.2.3.

**A3.3 No ambiguity [P mod FP Thm 4].** c_t' c_{-t} = Sum_{n>=1} (-1)^{n-1} p_n t^{n-1} for c_t = Sum lambda_i t^i, p_n = n! ch_n. Taking t^{2g-2} gives (2g-1)! ch_{2g-1}(𝔼) = (-1)^{g-1} lambda_g lambda_{g-1}. Further lambda_{g-1}^2 = 2 lambda_g lambda_{g-2}, so Int_{Mbar_g} lambda_g lambda_{g-1} lambda_{g-2} = (1/2) Int lambda_{g-1}^3 != 0 for g >= 2. So ch_{2k-1}(𝔼) != 0 on Mbar_k, k >= 2. For k = 1, Int_{Mbar_{1,1}} lambda_1 = 1/24.

**Colour isolation.** Work on Mbar_{k,1} with the insertion e_i. By the flat unit, (s_i.1)_{k,1} = p^*(s_i.1)_{k,0} and likewise for Hodge. p^*: H^*(Mbar_k;Q) -> H^*(Mbar_{k,1};Q) is injective for k >= 2 because p_*(psi_1 p^* x) = (2k-2) x, so ch_{2k-1}(𝔼) != 0 on Mbar_{k,1}. For k = 1, Mbar_{1,1} is used directly.

Now suppose k is minimal with a_k^{(i)} != b_k(wt_i). By A3.2 and A3.1, the degree 2k - 1 part of s_i.1 = Hodge_{wt_i} on Mbar_{k,1} reads (a_k^{(i)} - b_k(wt_i)) ((2k)!/B_{2k}) ch_{2k-1}(𝔼) = 0, a contradiction. So a_k^{(i)} = b_k(wt_i), that is s_i = h_{wt_i} = e^{-c_i}. This is BS (Normalisation), with BS's sign, over five weights.

**A3.4 Tensor structure [P].** For stable (g,n), Omega^Z = (eps1 eps2)^{g-1} Hodge_{(-eps1,-eps2)} Omega^{K_P2}. Hodge_v is multiplicative along boundary strata, including loops, because c_{-1/v}(𝔼) is (the normalised class, not Lambda^vee itself), and it is pulled back along p. So Hodge_v (R.omega) = R.(Hodge_v omega) = (R h_v).omega. Hence R^Z = R^{K_P2} h_{(-eps1,-eps2)}.

**A3.5 Oscillating integral [P mod draft SP1-SP5, Lemmas 2.1, 2.3-2.5].**

*Dependency.* A3.5 is not load-bearing for the graph sum: A1, A2, A3.1-A3.3 and A3.6 suffice for that. It is load-bearing only for Stage B's identity P_0 = Theta(z/eps1) Theta(z/eps2) Phi (B2). Its status there is P modulo SP1-SP5 and Lemmas 2.1, 2.3-2.5 of the draft.

*Field.* Asym has the constant c, c^2 = -1, and factors a^{-1/2} and sqrt(eps1 eps2), none of which need lie in K~w. The argument below runs over K^ := K~w(sqrt(-1), sqrt(eps1 eps2), a^{1/2} : a in {3 alpha_i, alpha_j - alpha_i, alpha_k - alpha_i}). Lemma 2.5(ii) goes through there unchanged: the theta-constants are the coefficient field, and Y, W are units, since Y_i(0) = -alpha_i and W_i(0) != 0 (the algebraic form of Lemma 2.5, stated after its proof in the draft). The colour-dependent branch sign is fixed together with the K~w-choice of (W_i^0)^{1/2}.

*Transfer of Lemma 2.4 to the formal ring.* Lemma 2.4 is stated pointwise on an analytic domain with y != 0, while the argument takes place in K^[[y]][[z]]. By Lemma 2.3, 𝒜 = c_0 W^{-1/2} Phi(alpha,Y;z), with c_0^2 = -1/3 and [z^n] Phi in Q[alpha,Y,W^{-1}]. Let 𝒫 := e^{-u/z} L_pf e^{u/z} with theta u = Y. The theta in 𝒫 acts on Q[alpha,Y,W^{-1}] through theta Y = -Y prod_j (Y+alpha_j)/W. So the z-coefficients of W^{1/2} 𝒫 (W^{-1/2} Phi) lie in Q[alpha,Y^{+-1},W^{-1}], because y = -prod_j (Y+alpha_j)/(27 Y^3) (Lhat is linear in y with coefficient 27 Y^3), and the fraction field of 𝒞 is Q(alpha,eps,Y). By Lemma 2.4 those coefficients vanish on an open set of (alpha,Y), hence they vanish identically as rational functions. Applying iota_i (Y -> Y_i in Kw[[y]], with Y_i and W(Y_i) units, since Y_i(0) = -alpha_i and W(-alpha_i) != 0, commuting with theta) gives 𝒫_i ((3W(Y_i))^{-1/2} iota_i(Phi)) = 0 in K^[[y]][[z]]. Here (3W(Y_i))^{1/2} lies in K^[[y]], because 3W(-alpha_i) is the product of the three a in the definition of K^. sqrt 3 itself is not adjoined.

*Kernel argument.* Put A_i := c (3W(Y_i))^{-1/2} iota_i(Phi), which equals c_0 W(Y_i)^{-1/2} iota_i(Phi) with c_0 = c/sqrt 3. Lemma 2.3 gives z-coefficients in K^[[y]], since W_i(0) != 0. By the above and A2.5, A_i and R~_{0i} are both formal solutions of 𝒫_i f = 0 with non-zero z^0 term. So by the kernel argument of Lemma 2.5(ii), A_i = kappa_i(z) R~_{0i}.

*Value at y = 0.* Evaluate Phi at Y = -alpha_i, an algebraic substitution with no limit, since W(-alpha_i) = alpha_i (alpha_j - alpha_i)(alpha_k - alpha_i) != 0. There H = diag(0,a,b) + (alpha_i/3)J and T_m = a xi_j^m/m! + b xi_k^m/m! + 3 alpha_i s^m/m!, with a = alpha_j - alpha_i, b = alpha_k - alpha_i, s = (xi_i + xi_j + xi_k)/3. In the linear coordinates (s, xi_j, xi_k) the phase splits, so SP1 (linear case) and SP4 factorise Phi into three one-dimensional expansions, with phases a(e^xi - xi), b(e^xi - xi), 3 alpha_i (e^s - s). With q = a e^xi (SP1 in dimension one), each factor is SP5 with beta = 0 and z -> -z, equivalently Stirling for Gamma(a/z): (-a)^{-1/2} Theta(-z/a). Hence A_i(0) = c prod_{a in {3 alpha_i, alpha_j - alpha_i, alpha_k - alpha_i}} a^{-1/2} Theta(-z/a), c = (-1)^{-3/2}.

*Conclusion.* With R~_{0i}(0) = (W_i^0)^{-1/2}, R_{0i} = R~_{0i} h_{wt_i} and Theta(z/eps) Theta(-z/eps) = 1, this yields A_i = c sqrt(eps1 eps2) Theta(-z/eps1) Theta(-z/eps2) R_{0i}. So eq:normNLR is derived and CCIT is not used. The rank-3 mechanism: h_{wt_i} divided by the Gamma limit is h_{(-eps1,-eps2)}, the same at all three colours. Unlike the conifold's nw, this is not 1 on T_F. On T_F the extra Theta(-z/E) comes from the SP5 factor at alpha = -E (B4), not from an h_nw cancellation.

**A3.6 Graph sum [P].** Combine A1.8, A2.7, A3.3 and the factorisation of Delta^{(2g_v-2+|H|)/2} as one Delta^{1/2} per half-edge and per dilaton leaf, with Delta^{g_v-1} left at the vertex. P_a^{(i)} is rational over Kw[[y]], because sqrt(W_i^0) cancels between Delta_i^{1/2} and R~. The Hodge-vertex form is R = R~ diag(h_{wt_i}), PPZ Prop 2.14, and diag(h).omega(e_i^n) = Delta_i^{g-1} Hodge_{wt_i}.

**Rank-3 specifics.**
- Three idempotents: at generic alpha all three Y_i are power series in y (simple roots of Lhat at y = 0), so nothing here involves t = y^{1/3}.
- Cube roots and conjugate colours arise only at T_F, where Lhat = (Y + E/3)^3 + 27 y Y^3 has a triple root at y = 0, Y_i = -(E/3)/(1 + 3 omega^i t), and Q = 0 is non-semisimple. That is why the constants must be fixed at generic alpha.
- O(-3): weights in A1.2. Omega^Z_{g,n} is valid for all g >= 0.

**Discharged for the draft.** (C2): (Normalisation) is now a theorem, A3.3, modulo the inputs of A2.8, A2.9, A3.2, A3.3. (C1) reduces to eq:J=I plus the QDE (A2.1-A2.4). (Symplecticity) is [P] (A2.6). Psi_{0i} = Delta_i^{-1/2} and Delta_i = 3 eps1 eps2 W_i are [P] (A2.5). Prop 2.6 is therefore P mod eq:J=I (CCIT15 Thm 31, s.1.2), the genus-0 QDE (derived in s.1.2) and SP1-SP5. BS 4056 is [P mod GP] (A1.3, B5).

---

## 3. Stage B: proof of Theorem FB_B

Theorem FB_A is used as stated. No finite-order computation is a proof step.

**B1. The ring 𝒞 [P].**
- Lhat has leading coefficient Delta_c, a unit of O, so 𝒞 is free on 1, Y, Y^2.
- On Lhat = 0: Y Lhat'(Y) = Y e_2(q) - 3 prod q = -W, with q_j = Y + alpha_j (direct expansion).
- Mod 𝔮: Lbar = (Y+E/3)^3 + 27 y Y^3. With t := -(Y+E/3)/(3Y) one gets t^3 = y and Y = -(E/3)/(1+3t).
- X^3 - y is irreducible over K_y (Eisenstein at y), so K_y[Y]/Lbar = K_y(t) = K(t) = L is a field. Thus 𝔮𝒞 is the unique maximal ideal.
- Wbar = E (Y+E/3)^2 != 0 in L, and Ybar != 0, Ybar+E/3 != 0, so W, Y, Y+alpha_j are units of the local ring 𝒞 (their residues are non-zero).
- theta = y d/dy preserves O, because theta 𝔮 lies in 𝔮 and denominators stay outside 𝔮. It extends to 𝒞 by theta Y = -27 y Y^3/Lhat'(Y) = -Y prod q_j/W. It reduces to the unique extension of y d/dy to L, namely (t/3) d/dt. Check by hand: -Y(Y+E/3)/E = (E/3) t/(1+3t)^2 = (t/3) dY/dt.
- Lhat is linear in y with coefficient 27 Y^3, so 𝒞 is a domain with fraction field Q(alpha,eps,Y). This is the fact used in the transfer of A3.5.

*Colours.* At generic alpha, Lhat has three simple roots Y_j in Kw[[y]], Y_j(0) = -alpha_j (Hensel; Stage A's Y_j). Non-zero polynomials in y are units of Kw((y)), so O embeds in Kw((y)), and iota_j: 𝒞 -> Kw((y)) sends Y -> Y_j, commuting with theta. (iota_0,iota_1,iota_2) is an isomorphism after (x)_O Kw((y)), so Tr_{𝒞/O} = Sum_j iota_j. For T in 𝒞^{(x)n}, Tr^{(x)n}(T) = Sum_c (iota_{c_1} (x) ... (x) iota_{c_n})(T) lies in O.

*Reduction.* The trace of a free module commutes with base change, so the reduction of Tr^{(x)n}(T) is Tr_{L/K_y}^{(x)n}(pi^{(x)n} T). After adjoining omega this equals Sum_c of pi^{(x)n} T evaluated at t_v = omega^{c_v} t. So the ring is 𝒞 (a local etale O-algebra, because Lhat'(Y) = -W/Y is a unit), no roots of the cubic are needed, and Galois invariance is replaced by the trace. At generic alpha 𝒞 is not Galois over O.

**B2. Row 0 at generic alpha [P, given A3.5 and draft Lemma 2.3].**
- Stage A, in the form written in A3.5: c (3W(Y_i))^{-1/2} iota_i(Phi) = c N R_{0i} in K^[[y]][[z]], where c (3W(Y_i))^{-1/2} = c_0 W(Y_i)^{-1/2} with c_0 = c/sqrt 3. The transfer step of A3.5 is what makes the oscillating-integral statement an identity in this formal ring.
- By formal stationary phase in log coordinates (eq:spotP2), Asym_i = c_0 W_i^{-1/2} Phi(alpha,Y_i;z), with c_0 a universal constant. The jets of the phase at the critical point are q_j delta - Y/3^{m-1}, polynomial in (alpha,Y), and H^{-1} = (3/W) adj H.
- The z^n coefficient is a finite sum over k <= 2n factors T_m (draft Lemma 2.3), and it depends on y only through Y. This remains valid at y = 0, where q_i = 0, because W_i(0) != 0.
- With Delta_i = 3 eps1 eps2 W_i and N = sqrt(eps1 eps2) Theta(-z/eps1) Theta(-z/eps2): Delta_i^{1/2} R_{0i} = sqrt(3)(c_0/c) Theta(z/eps1) Theta(z/eps2) Phi.
- The z^0 term (Delta^{1/2} Psi_{0i} = 1) gives sqrt(3) c_0/c = 1. All branch choices cancel, so P_0^{(i)} = iota_i(P_0), with [z^n] P_0 in 𝒞.
- Note that h_{wt_i} and R~ are each singular on the diagonal, while the product is not: h_{wt_i} = h_{(3 alpha_i, alpha_j - alpha_i, alpha_k - alpha_i)} h_{(-eps1,-eps2)}, and the Gamma part is absorbed in Phi.

**B3. Rows 1, 2 and the pairing [P, given A2.6, eq:RHb as derived in Stage A, draft Prop 2.12].**
- From R_1 = (z theta + Y) R_0/A, R_2 = (z theta + Y - e_1 E(y)) R_1/D, and P_a := Delta^{1/2} R_a: P_1 = 𝐐/A, P_2 = D^{-1} A^{-1} [𝐓𝐐 - z F 𝐐 - e_1 E(y) 𝐐].
- By Prop 2.12: 1/(AD) = Delta_c A and e_1 E(y) = (e_1 X/2)(A^{-2} - 1).
- With eta^{11} = -3 e_1 eps1 eps2, eta^{12} = -3 eps1 eps2, eta^{22} = 0 (from eta_00 = 1/(3 e_3 eps1 eps2), eta_12 = -1/(3 eps1 eps2), eta_22 = e_1/(3 eps1 eps2), which also follow from the residue formula), the A^{-2} terms cancel between eta^{11} P_1 P_1' and eta^{12}(P_1 P_2' + P_2 P_1'). This gives (Pair). It is the frame-free reason for BS prop:fingenTLE, and it agrees with BS eq:dFR.
- Edges: N_ij(psi,w) := <P^{(i)}(-psi), P^{(j)}(-w)> satisfies N_ij(-w,w) = delta_ij Delta_i (symplecticity, A2.6). Hence [delta_ij Delta_i - N_ij(psi,w)]/(psi+w) = -Sum_{M>=1,b} n_{M,b} w^b (psi^M - (-w)^M)/(psi+w), and the M = 0 terms, which contain delta_ij, never enter. So E_ij(a,b) = -Sum_{k=0}^b (-1)^k n_{a+k+1,b-k}, a tensor in 𝒞 (x)_O 𝒞[𝖥], with no idempotent needed. Fix an orientation of every edge to define T_{Gamma,dec}. The map 𝒞 (x)_O 𝒞 -> Kw((y))^9 given by the iota_i (x) iota_j is injective, since 𝒞 (x)_O 𝒞 is free over O with basis Y^a (x) Y^b, 0 <= a,b <= 2, O embeds in Kw((y)), and the matrix (Y_i^a Y_j^b) is the Kronecker square of a Vandermonde matrix, invertible because the Y_i are distinct. Hence the tensor of E(a,b) equals the flip of the tensor of E(b,a) in 𝒞 (x)_O 𝒞[𝖥], T_{Gamma,dec} does not depend on the orientation, and this passes to L (x)_{K_y} L under pi. For a loop, apply the multiplication 𝒞 (x)_O 𝒞 -> 𝒞.
- **𝖥 free.** The telescoped E_ij(a,b) is defined for 𝖥 an indeterminate, and the delta form holds there too. The 𝖥-part of Sum_{M+b=n} (-1)^b [z^M w^b]<P(z),P(w)> is

      3 eps1 eps2 Delta_c 𝖥 Sum_{i+j=n-1} [(-1)^j + (-1)^{j+1}] 𝐐_i (x) 𝐐_j = 0.

  Equivalently, the 𝖥 term of <P(z),P(w)> carries the factor (z+w). So symplecticity and the delta form hold with 𝖥 free. No theta F occurs in leg-free terms, so Stage C may treat 𝖥 as a free symbol.
- Delta^{g_v-1} is in 𝒞, including g_v = 0 (Delta is a unit). Dil(d) = (-1)^d [z^{d-1}] P_0.
- The decoration sum is finite (Sum a + Sum d = dim, d >= 2). This gives (a) and (c) at generic alpha, with |E| <= 3g-3 bounding deg_𝖥.

**B4. Reduction of the rows [P, given draft Theorem A, Cor 2.2, Lemma 2.4, SP1-SP5].**

pi(Phi) is a power series in z with coefficients in K[Y, Wbar^{-1}]. It is computed in the s-chart (draft Thm 2.7(e), steps 1-3); the derivation is written here, because the draft does not show that the s_3-dependence is a pure power.

*Chart.* At T_F, in the coordinates a_j := q_j^{1/3} and s_n := e_n(a), the critical point of colour 0 is a = c'(1, omega, omega^2), with distinct entries, so s is a local biholomorphism there (SP1). With kappa := 3 + 1/t = E v, the phase is

    s_1^3 - 3 s_1 s_2 + kappa s_3 - E log s_3,

which splits exactly, with critical point s_1 = s_2 = 0, s_3 = 1/v and block-diagonal Hessian (SP4). The measure is dOmega = (1/3) d^3x = 9 ds/(s_3 V), with V^2 = Disc, and

    Disc/(-27 s_3^2) = 1 - (2/3) s_1 s_2/s_3 + (4/27) s_1^3/s_3 + (4/27) s_2^3/s_3^2 - (1/27) s_1^2 s_2^2/s_3^2,

so only integer powers of 1/s_3 occur, and 1/V contributes one factor of 1/s_3 times a series in these monomials.

*The (s_1,s_2)-block.* The only non-zero contraction is <s_1 s_2> = z/3. All monomials in the expansion have weight 0 for wt(s_n) = n. So a term s_1^a s_2^b s_3^{-n} of the amplitude, with k insertions of s_1^3/z, has a + 2b = 3n and needs a + 3k = b. It carries z^{b-k} and n = b - k. Put l := b - k. The total of the (s_1,s_2)-block is Sum_l r^W_l z^l s_3^{-l-2} ds_3, up to the constant Gaussian factor. This defines BS's Wick numbers r^W_l, with r^W_0 = 1 after normalisation.

*The s_3-integral.* Substitute s_3 = -q/kappa. The phase becomes -(q + E log q) up to a constant, and s_3^{-l-2} ds_3 = (-kappa)^{l+1} q^{-l-2} dq. So SP5 applies with alpha = -E, beta = -l-1, giving (-E)^{-l-3/2} Theta(-z/E, -l-1), where Theta(w, -m) = Theta(w)/prod_{j<=m}(1 - jw). The factors combine as (-kappa)^{l+1} (-E)^{-l-1} = v^{l+1}, because kappa/E = v, with no l-dependent sign. So the l-th term is r^W_l z^l v^{l+1} Theta(-z/E)/prod_{j=1}^{l+1}(1 + jz/E), up to a constant, and

    𝒜 = K_0 Theta(-z/E) v S_Wick(v,z),

with K_0 independent of z and y.

*Recursion.* By Lemma 2.4, e^{u/z} 𝒜 solves eq:pf. Theta(-z/E) is invertible and z-only, so 𝒫(v S_Wick) = 0. Corollary 2.2 with r^W_0 = 1 gives r^W_l = r_l for all l.

*Conclusion.* With Wbar^{-1/2} proportional to v and Phi_0 = 1: pi(Phi) = Theta(-z/E) S(v,z). This is an identity in the single rational variable Y (v = 3/(E+3Y)), so it holds at all three roots without a chart or monodromy discussion; this replaces step 5 of the draft's proof. Hence pi(P_0) = Theta(z/eps1) Theta(z/eps2) Theta(-z/E) S = h_w S, and "the z-only factor is 1" is now a consequence: the constants enter only through Stage A's A3.5. pi(𝐓) = 𝐓bar, since pi(W) = E (Y+E/3)^2 and pi(Delta) = Delta_i; pi(eta^{ab}) = etabar^{ab}. By B1, pi commutes with theta. Thus pi^{(x)V}(T_{Gamma,dec}) is the T_F coloured term. Its telescoped edge equals the delta-form edge by draft Prop 2.10 (proved from the closed form alone) and B3.

*Against draft Thm 2.7.* Thm 2.7(a)-(c) is correct and gives regularity locally near (y0,alpha*). But it is stated for germs continued along paths from y near 0 (Prop 2.6), and BS's 4196-4198 ansatz is not needed. What the identification with Stage A's series needs is B2: rationality in (alpha,Y) with denominator a power of W.

*Analytic corollary [P].* For y0 not in {0,-1/27} and (y,alpha) near (y0,alpha*), the root Y near Y_i is analytic (Lhat' = -W/Y != 0), and all factors are rational in (alpha,eps,y,Y) with no pole there, at every z-order. This is the all-orders replacement for "first few n" and "by continuity" (BS 4347, 4198).

**B5. GW side [P mod Graber-Pandharipande; A1.2].**

*Graber-Pandharipande for T.* GP's formula (1) and Appendices A, C are stated for a C^*-action. Apply them to a generic one-parameter subgroup through which T acts with the same fixed loci. Both sides are rational in (alpha,eps), and such subgroups are Zariski dense, so the formula holds for T.

*d >= 1.* Let M = Mbar_g(P2,d), proper, with T = (C*)^3 x (C*)^2 acting through alpha (a linearisation with H|P_i = -alpha_i, O(-3)|P_i = 3 alpha_i). H^0(C, f^*O(-3)) = 0: a section vanishes on non-contracted components, is constant on contracted ones, and C is connected. So R^1 := R^1 pi_* f^* O(-3)_{3 alpha} is a vector bundle of rank 3d+g-1, and e_T(R pi_* f^* N)^{-1} = e_T(R^1) Lambda^vee(-eps1) Lambda^vee(-eps2)/(eps1 eps2), an honest class on M after inverting eps1 eps2. By A1.2 (all weights of R pi_* f^* N move on fixed loci, and [F]^vir and N^vir compare) and GP on M, BS's localisation definition equals Int_{[M]^vir_T} of this class. The T-equivariant pushforward lands in Q[alpha,eps], which gives the polynomiality. The degree is (3d+g-1) + 2g - (3d+g-1) = 2g, and the expression is symmetric in eps1, eps2.

*d = 0, g >= 2.* The same argument on Mbar_g x P2, with e_T(O(-3)_{3 alpha}) = -3H inverted as in the statement, using H^3 + e_1 H^2 + e_2 H + e_3 = 0.

*Passage to T_F.* Restriction to T_F (on a finite cover; rational coefficients make this harmless) is alpha_i -> E/3 and is compatible with pushforward. On T_F: Mbar(Z,d)^{T_F} = M (a compact curve maps to a point of the affinisation, so for T_F-fixed maps it lies in the zero section), fixed part = obstruction theory of P2, N^vir = R pi_* f^* N with weights E, -eps1, -eps2 non-zero. At alpha*: H = H' - E/3, and (H^2+e_1 H+e_2)/(3 e_3) = 1/(E - 3H') mod H'^3 (checked by hand). e_3(alpha*) = (E/3)^3 is not in the prime. This proves BS 4056 and (Contravariance) here.

*Constant term.* Here H' is the non-equivariant hyperplane class. Mod H'^3, c_top(𝔼^vee box T_P2) = -3H' lambda_g lambda_{g-1} + 9 H'^2 lambda_g lambda_{g-2}, using x1+x2 = 3H', x1 x2 = 3H'^2, lambda_g^2 = 0, lambda_{g-1}^2 = 2 lambda_g lambda_{g-2} (from c(𝔼) c(𝔼^vee) = 1 [C]). In P(H') the term with 1/(E-3H') carries lambda_g, so it dies against lambda_g in the prefactor: eps1 eps2 GW^{T_F}_{g,0} is a polynomial, with no 1/E. For g = 2 this gives 9 E^2 Int lambda_1 lambda_2 = E^2/640. For g = 3, put u = E - 3H', p = eps1 eps2. Mod lambda_3 and with lambda_1^2 = 2 lambda_2, P(H') = (u^2 - lambda_1 u + lambda_2)(p^2 + lambda_1 p E + lambda_2 E^2 + ...), so [P(0)]_2 = lambda_2 (E^2 - p)^2 and [P'(0)]_1 = -3 lambda_1 (2 E^2 p - p^2). Hence GW^{T_F}_{3,0} = 9 I [(E^2 - p)^2 + 2 E^2 p - p^2] = 9 I E^4, with I = Int lambda_1 lambda_2 lambda_3 = 1/1451520 (FP Theorem 4).

**B6. The num/den step [P].**
- Write 𝒯_g = Sum_k (p_k/den) 𝖥^k, with p_k, den in 𝔸 and den not in 𝔮.
- By Stage A and B3, in Kw[[y]]: den F^T_g = Sum_k p_k F^k.
- Every term lies in the subring O_w[[y]]: F^T_g = Sum Q(y)^d GW^T_{g,d} with Q(y) in y Q[[y]] (mirror map, independent of alpha) and GW^T in O_w by B5; F in Q[[y]]; p_k, den polynomial.
- Apply the ring map rho: rho(den) F_g^{T_F} = Sum rho(p_k) F^k in the domain K[[y]].
- rho(den) = den|_{alpha*} != 0 in K[y], since den is not in 𝔮. Dividing in K((y)) gives F_g^{T_F} = (reduction of 𝒯_g)(F).
- By B1 and B4 the reduction of 𝒯_g is the T_F colour sum, also coefficientwise in the indeterminate 𝖥, so the reduction of the polynomial is Tbar_g[𝖥] of s.4. This is (d).
- rho(den)(0) may vanish, and this does no harm. Model case: Norm_{𝒞/O}(W) = prod_i W(Y_i) = -e_3 Delta_c^2 prod_{i<j}(Y_i - Y_j)^2. At y = 0 it is -e_3 prod_{i<j}(alpha_i - alpha_j)^2, while on T_F it equals E^9 y^2/Delta_c^2 (by hand: prod_i (Y_i+E/3) = y E^3/Delta_c), which vanishes at y = 0.

Why no limit is taken. Coefficientwise in y, a colourwise factor blows up on the diagonal. At fixed y0 the y-series of Y_j(y) diverges for alpha near alpha*, because the branch points of Y(y) tend to y = 0. Hence neither a pointwise limit of the y-series nor a colourwise limit is valid. The identity is proved by a ring map, with no topology.

Draft Prop 2.11's phi_{k,m}, P_{k,m}/Q arose from a posteriori finite generation, and non-divisibility by delta was deduced by contradiction. The round-5 critic item S1(c) is right that the deduction needs independence of A, F over functions of y: at fixed y0 the graph sum is one number, so a pole of some phi need not show. Here den not in 𝔮 is built in, because the coefficient of 𝖥^k is a multi-trace of an element of 𝒞^{(x)n}. So no independence is used, only one variable 𝖥 is present (A has cancelled), and no rational torus is needed. The square-root reason for rational tori disappears, since normalised rows have no square roots. eq:adelta is therefore unnecessary.

**Proof of (e) [P].**

*Coefficient ring.* [z^n] of Pbar_0, Qbar, 𝐓bar Qbar are Laurent polynomials in v, hence in K[t^{+-1},(1+3t)^{+-1}]. The colour sums lie in K(omega)[t^{+-1},(1+3 omega^i t)^{-1}], which is integral over R_omega := K(omega)[y^{+-1},Delta_c^{-1}], because t^{-1} = t^2/y and (1+3 omega^i t)^{-1} = prod_{k != i}(1+3 omega^k t)/Delta_c. The colour sum lies in K_y (B1). R_omega is integrally closed, so the sum lies in R_omega ∩ K_y. K(omega)/K is a constant-field extension, so denominators are detected by the y-adic and (1+27y)-adic valuations over K[y] (a gcd computed over K stays a gcd over K(omega)). So R_omega ∩ K_y = K[y^{+-1},Delta_c^{-1}].

*F is transcendental over ℂ(y).* No modularity is used. A solves L2b = Delta_c theta^2 + 27 y theta + 6y, with exponents (0,0) at Delta_c = 0. Let wc be the analytic solution with wc(0) = 1. The local solutions are spanned by wc and wc log xi + g_1, with g_1 analytic. Since a_d 27^{-d} ~ sqrt3/(2 pi d) > 0, A -> +infinity as y -> -1/27 from the right. So A = a wc + b (wc log xi + g_1), with b != 0. (Explicitly b = -sqrt3/(2 pi), from the formula of C Step 6 with u -> 1-u.) Continuing n times around xi = 0 gives A_n = A + 2 pi i n b wc. If F were algebraic over ℂ(y), it would have finitely many branches, so theta log(A_n/A_m) = 0 for some n != m, that is A_n = c A_m with c constant. For c = 1 this contradicts b wc != 0. For c != 1 it makes A proportional to wc, contradicting b != 0. Hence F is transcendental over ℂ(y).

So F is transcendental over ℂ(t), since t is algebraic over ℂ(y), and, taking coefficients in eps as in C Step 4(i), over K_y. So the polynomial in 𝖥 is unique.

BFGW Lem 4.2 (equivalently prop:algindABC with eq:XFABC and eq:ABCmod) gives more, the independence of A and F, and is needed only for statements involving A (C Step 4).

*A consequence for Stage C.* The rows of draft Def 1.3 depend on c1 only as (Pair) with 𝖥 -> XX (draft Prop 2.10 S0, for the same P_2 decomposition), so Cf acts graph by graph as 𝖥 -> F_cf.

**Rank-3 specifics, as handled.** Cube roots t = y^{1/3}: only in the residue field L. Complex-conjugate regular colours: in (d) they are t_1, t_2, conjugate when t and eps are real, with no special treatment, since the sum is a trace. Three idempotents at generic alpha: iota_0, iota_1, iota_2. Non-semisimplicity at Q = 0 on T_F: the constants are fixed at generic alpha (A3) and then carried by the ring map; Teleman is not applied on T_F. Triple root at y = 0: visible only as t-poles of coloured terms, which cancel in the total as a corollary of (d) (the left side is in K[[y]]). Regularity at y = 0 is not an input.

---

## 4. Stage C: proof of the frame change, Proposition 13.2_P2, Lemma (t_cf), Corollaries G and G_{k,h}

Stages A and B are used as stated. No finite-order computation is a proof step.

**Notation.** For a colour i, the rows of draft Def 1.3 in normalised form, with T_i := z theta + Y_i - (z/2) theta log Delta_i, are: Phat_0^{(i)} := S(v_i,z); Q_i := T_i S(v_i,z); Phat_1 := Q_i/c1; Phat_2 := Delta_c c1 [T_i Q_i - z 𝖥 Q_i - d2 Q_i], with d2 := E(c1^{-2}-1)/(2 Delta_c). Here c1 and 𝖥 are independent symbols. G_g[𝖥] is the decoration sum (0.1) of the draft for these rows, with vertex classes Delta_{c_v}^{g_v-1} Omega_w. Tbar_g[𝖥] is the same sum with rows h_w(z) Phat_a and point vertices (Stage B's 𝒯bar_g, with 𝖥 free).

**C Step 0. The abstract setting [P].** Let L' := L_omega(kappa_0)(c1,𝖥), kappa_0^2 = 3 eps1 eps2 E. It is finitely generated over Q, hence embeds in ℂ, so PPZ's statements apply and identities pull back. Put V := L'^3 with basis H^0 = 1, H^1, H^2 and eta^{ab} of draft s.1.3. Normalise Delta_i^{1/2} := kappa_0/v_i.

(Pair) [P, algebra, no theta on c1]. Expanding with eta^11 = E eta^12, eta^22 = 0, and 2 Delta_c d2 = E(c1^{-2} - 1), the c1^{-2} terms cancel and

    <Phat^i(z),Phat^j(w)> = eta^00 S_i S_j' - 3 eps1 eps2 E Q_i Q_j' - 3 eps1 eps2 Delta_c (Q_i (T Q_j)' + (T Q_i) Q_j') + 3 eps1 eps2 Delta_c (z+w) 𝖥 Q_i Q_j',

where a prime means evaluation at w. This agrees with Stage B (Pair), and it is an exact identity in L'[[z,w]]. At w = -z the term in 𝖥 and all dependence on c1 disappear. So S0 and S1-S6 of draft Prop 2.10 give <Phat^i(z),Phat^j(-z)> = delta_ij Delta_i in L'. The proof of Prop 2.10 uses c1 only through S0, so it holds with c1 and 𝖥 free; B3 gives the same statement for the 𝖥-part.

Hence Psi_{ai} := Delta_i^{-1/2} Phat_a(0) satisfies Psi^T eta^{-1} Psi = 1. Define f_i by eta(H^a,f_i) = Psi_{ai}. Then eta(f_i,f_j) = delta_ij and, since Phat_0(0) = 1, 1 = Sum_i Delta_i^{-1/2} f_i. Define R'(z) by eta(H^a,R'(z) f_i) = Delta_i^{-1/2} Phat_a^{(i)}(z). Then R'(0) = 1 and R' is symplectic. omega(e_i^n) := Delta_i^{g-1}, with e_i = Delta_i^{-1/2} f_i, is a semisimple TFT with unit 1. From symplecticity, eta(f_i,R'^{-1}(z) H^a) = eta(R'(-z) f_i,H^a) = R'_{ai}(-z).

**C Step 1. Unfolding lemma (C2) [P, PPZ Defs 2.2, 2.5, 2.13].** Let Theta_0 be a rank-1 CohFT with unit, pulled back under forgetful maps, and let R be symplectic with rows r_{ai}.
- T(z) = z(1 - R^{-1}(z)1) = z Sum_i (Delta_i^{-1/2} - r_{0i}(-z)) f_i = Sum_i e_i z(1 - Delta_i^{1/2} r_{0i}(-z)). This is O(z^2), because R(0) = 1.
- The coefficient of e_i z^d is (-1)^d [z^{d-1}] Delta_i^{1/2} r_{0i}. For d = 1 it vanishes in T. The term -z of q_t is the -psi of the dilaton shift, so q_t = T^{(i)} - psi exactly, where T^{(i)} is the coefficient of e_i in T(psi), not the operator T_i of the Notation.
- The bivector coefficient on f_i (x) f_j is [delta_ij - Sum r_{ai}(-psi) eta^{ab} r_{bj}(-w)]/(psi+w). On e_i (x) e_j, multiply by Delta_i^{1/2} Delta_j^{1/2}, which gives (1.1) of the draft.
- (Theta_0 (x) omega)(e_i^n) = Delta_i^{g-1} Theta_0, and mixed insertions vanish.

Therefore Int_{Mbar_g} (R.(Theta_0 (x) omega))_{g,0} equals the sum (0.1) with Theta_0 at the vertices, m_v-tuples weighted by 1/m_v!, and d >= 2.

*Frame-B leaves.* The leaf series proper is T'_i = psi - psi S(v_i,-psi), starting at psi^2. The series q starts at psi^1. At colour 0, t_c = 0: S(0,-psi) = E/(E-psi), so T' = -psi^2/(E-psi) = q_0 + psi, exactly the T of draft s.1.5. The rescaled rho q_t has linear coefficient -rho; it is handled by the domain 𝔇 and the dilaton equation, not by the decoration sum. These are consistent.

Rows with Theta_0 = 1, R = R'h_w give Tbar_g[𝖥]. Rows with Theta_0 = Omega_w, R = R' give G_g[𝖥].

**C Step 2. Frame change (C1) [P mod draft Prop 2.10, A3.2 (Mumford as in FP), PPZ Defs 2.2, 2.5, 2.13, Prop 2.14. Stages A, B are needed only to identify Tbar_g[𝖥] with the reduction of 𝒯_g, B6, not for this identity.]**
- A3.2: Hodge_v = h_v.1 with h_v = prod_a Theta(-z/a). For v = w = (E,-eps1,-eps2) this is h_w = Theta(-z/E) Theta(z/eps1) Theta(z/eps2). So the choice of w and the sign are not fitted.
- For the scalar h_w: graphs are monochromatic, Sum_v (g_v - 1) + |E| = g - 1, and the translation has coefficient 1 on e_i. So (h_w.omega)(e_i^n) = Delta_i^{g-1} (h_w.1)_{g,n}, with mixed insertions zero. That is, h_w.omega = Omega_w (x) omega as CohFTs with unit, including g = 0, where Omega_w = 1.
- PPZ Prop 2.14: (R' h_w).omega = R'.(h_w.omega) = R'.(Omega_w (x) omega). The dictionary gives (R' h_w)_{ai} = R'_{ai}(z) h_w(z), so large-radius rows are h_w(z) Phat_a, which are Stage B's Pbar with 𝖥.
- Integrate over Mbar_g and use C Step 1: Tbar_g[𝖥] = G_g[𝖥] in L'. Both sides lie in L_omega[𝖥]: by (Pair) of C Step 0 the edges do not contain c1, Dil uses only S, and by Remark (signs) of the draft kappa_0 does not survive.
- h_w contains only z and eps, and log h_w is odd in z. So h_w is t-independent, colour-independent, and h_w(z) h_w(-z) = 1. For the identification with the reduction of 𝒯_g, which is the only use of Stage B in this step: by B4, pi(P_0^{(i)}) = h_w S(v_i,z) as an identity in the single rational variable Y, hence with the same h_w at all three roots. So Stage B proves the point-frame rows are h_w Phat, not merely up to a colour-dependent z-only factor.
- No rank-3 input enters except Prop 2.10: the frame change is scalar. The round-5 critic item S1(a) is correct. The worry in draft s.8.3 about Lemma 13.3 is void, since h_{wt_i} is singular on the diagonal but Stage B's Phi absorbs it. On T_F the weights w are exactly the three non-zero normal weights of P2.

**C Step 3. Dictionary at c1 = A [P, draft Prop 2.12].** Delta_c A^2 = 1/D and d2 = e_1(alpha*) E(y) make Def 1.3 at c1 = A exactly eq:RHb on rows 1 and 2 (theta c1^{-1} = -c1^{-1} XX, so Phat_2 = Delta_c c1 [T T S - z XX T S - d2 T S]). This equals Stage B's Pbar_2 after h_w. Delta^{1/2} D_{i,n} = [z^{n-1}](1 - Pbar_0(-z)) = Dilbar_i(n), and eq:edgeterm times Delta_i^{1/2} Delta_j^{1/2} is Ebar. Stage B (d) reads F_g^{T_F} = Tbar_g[F] in K((y)).

**C Step 4. Definition and well-definedness of Cf [P for the leg-free part; statements involving A are P mod BFGW Lem 4.2].**

(i) *Independence.* Only the transcendence of F is used on leg-free terms, because A cancels. F is transcendental over ℂ(y) by the proof of (e) in s.3 (no eq:ABCmod). Now choose a branch of t on U; t is algebraic over ℂ(y), so F is transcendental over ℂ(t). For eps, eps3: the target of ev_A is M(U)(eps1,eps2)(eps3), with basis {1, eps3} over M(U)(eps1,eps2), since E = eps1 + eps2 is not a square there. If a polynomial in F with coefficients in ℂ(t)(eps1,eps2,eps3) vanished, take components along {1, eps3}, clear denominators, and take the coefficient of each monomial in eps. This gives a relation over ℂ(t). So F is transcendental over L~.

For statements involving A (ev_A on all of 𝐆, Cf(R_{ai}), the legged form of eq:gtcf), cite BFGW Lem 4.2: S, X, I11 are algebraically independent over ℂ, with I11 = A, S = F - (X-1)/3. So A, X, F are independent, and no birational step through BS's B, C is needed. tr.deg_ℂ ℂ(y,A,F) = 3 on U, and the same argument gives A, F algebraically independent over L~.

(ii) *Ring.* BS's G_{T_F} has abstract generators Y_i with relations L(Y_i), so coloured terms do lie in it abstractly, as polynomials in 1/(3Y_i+E), Y_i^{+-1}. It is not a domain: it has zero-divisor components (for example factors of Y_0 - Y_1). The specialisation to functions does not land in Q_T((y)) on T_F, because the roots need t. iota is the component with the labelled roots; any other labelling differs by S_3, and the invariant total is unaffected. Hence:
- coloured Cf is ev_cf on the lift in 𝐆 with t a formal variable, and this is how it is used in C Step 5;
- on functions, Cf := ev_cf o ev_A^{-1} on ev_A(𝐆). For a single coloured function this depends on the branch of t on U. On the colour-summed total it is independent of that branch and of lifts: a difference of lifts lies in ker(iota)[𝖠^{+-1},𝖥], which ev_cf also kills.

No analytic continuation of t from y = 0 to the conifold is used. ev_cf embeds L~ by the branch at xi = 0, which defines colour 0. This is the branch of draft s.0.1 and Lemma 1.2(b), with t_i = omega^i t and t -> -1/3. There v_0 = (1+3t)/(E t) = -xi/E + O(xi^2), v_1 and v_2 are non-zero at xi = 0, and Y_0 is the expression at BS 5340. So the colour called 0 under ev_cf is the singular colour of Theorems LC_P2, R_P2 and Corollary SV_P2. The identification is needed only to match per-colour statements, since the colour-summed total is independent of it. For real t, eps, colours 1 and 2 are conjugate and their terms are conjugate.

(iii) *theta.* theta A = A F and theta F = -F^2 - 27 y X F - 6 y X (from L2b A = 0). A_cf satisfies L2b, so F_cf satisfies the same Riccati equation. So Cf commutes with theta. Hence, with sqrt(3 eps1 eps2) adjoined to L~ (R_{ai} contain Delta_i^{-1/2} = (3 eps1 eps2)^{-1/2} W_i^{-1/2}), Cf(R_{ai}) are the rows of Def 1.3 at c1 = A_cf times h_w, and Cf(D) = X/A_cf^2, as in BS 5172. Up to the constant GW^{T_F}_{g,0} (BS's left side is printed with d > 0), this establishes BS eq:gtcf as a theorem.

(iv) *Uniqueness of the lift.* K[y^{+-1},Delta_c^{-1}] = K[X^{+-1},(1-X)^{-1}], the image of G^{S_3}_{T_F} (eq:GTFS3), using y = (1-X)/(27X). By (i) the polynomial in 𝖥 is unique, and finite generation is obtained from Stage B rather than cited. Transcendence is needed only to make Cf well defined on functions. The identities themselves are identities in L_omega[𝖥] by C Steps 0-2 and Stage B (d).

**C Step 5. Proposition 13.2_P2 [P mod the inputs in its label, s.0.2].** Cf(F_g^{T_F}) = ev_cf(Tbar_g[𝖥]) = ev_cf(G_g[𝖥]) = G_g[F_cf], by C Step 2. By C Step 0's (Pair), the edges of (0.1) at c1 = C1cf have 𝖥 = XX = theta log C1cf, c1 drops out, and Dil does not involve c1. With the logarithmic form in the Lemma below (C Step 6, no Barnes), A_cf = C1cf and F_cf = XX, so G_g[F_cf] is (0.1) at c1 = C1cf, term by term. The d = 0 term GW^{T_F}_{g,0} lies in K and is fixed by Cf.

**C Step 6. Lemma (t_cf) [P mod DLMF 15.6.7, Barnes's second lemma, Gauss].** Put u = -27y, so xi = 1-u, and w(u) := 2F1(1/3,2/3;1;1-u) = wc. (L2b is the hypergeometric equation with a = 1/3, b = 2/3, c = 1, in u; it is invariant under u -> 1-u, and A(u) = 2F1(1/3,2/3;1;u) = w(1-u).)

*Barnes representation.* DLMF 15.6.7 with c = 1 gives, for u in (0,1],

    w(u) = K_B (1/2pi i) Int Gamma(-s)^2 Gamma(1/3+s) Gamma(2/3+s) u^s ds,   K_B = (sqrt3/2pi)^2,

on a contour with -1/3 < Re s < 0, separating the poles of Gamma(1/3+s), Gamma(2/3+s) from those of Gamma(-s).

*Logarithmic form.* Closing right, the double poles at s = n and Gauss's psi(1/3)+psi(2/3)-2psi(1) = -3 log 3 give w = -(sqrt3/2pi)(A l + h), with l = log(-y), and h = Sum_d 3 a_d (H_{3d} - H_d)(-y)^d as in draft Prop 2.12 Step 0. This uses Sum_{k<n}(1/(k+1/3)+1/(k+2/3)) - 2H_n = 3(H_{3n} - H_n), and log u = l + 3 log 3.

*Derivative.* From eq:Ian: I_1^{[2]} = (log^2 y)/2 + g log y + g2 and theta I_1^{[2]} = A log y + h (checked by hand from the printed series), with I_0^{[1]} = log y + g, A = 1 + theta g, h = g + theta g2. With log y = l + i pi:

    I_1^{[2]} - i pi I_0^{[1]} - pi^2 = l^2/2 + g l + g2 - pi^2/2,

and so theta t_cf = (sqrt3/2pi)(A l + h) = -w. So A_cf = -w = C1cf and F_cf = XX.

*Constant.* With t_c as in draft s.1.2 (theta t_c = -wc, t_c(xi = 0) = 0), t_c = T(u) - T(1), where

    T(u) := K_B (1/2pi i) Int Gamma(-s)^2 Gamma(1/3+s) Gamma(2/3+s) (-u^s/s) ds,   -1/3 < Re s < 0,

so that s = 0 lies to the right of the contour. This contour is load-bearing: the triple pole at s = 0 must be the one crossed when closing right. The integrand decays like e^{-2 pi |Im s|}, up to powers, so T is continuous at u = 1 and theta_u may be taken under the integral; thus theta T = -w.
- *Value at u = 1.* Since -1/s = Gamma(-s)/Gamma(1-s), T(1) is K_B times the Barnes integral of Gamma(-s)^3 Gamma(1/3+s) Gamma(2/3+s)/Gamma(1-s). Barnes's second lemma, in the printed form

      (1/2pi i) Int Gamma(a'+sigma) Gamma(b'+sigma) Gamma(c'+sigma) Gamma(1-d'-sigma) Gamma(-sigma)/Gamma(e'+sigma) dsigma
         = Gamma(a') Gamma(b') Gamma(c') Gamma(1-d'+a') Gamma(1-d'+b') Gamma(1-d'+c')/(Gamma(e'-a') Gamma(e'-b') Gamma(e'-c')),  e' = a'+b'+c'-d'+1,

  applies with s = -(sigma + 2/3) and (a',b',c',d',e') = (2/3,2/3,2/3,4/3,5/3) (equivalently, the symmetric form with parameters (0,0,0,1/3,2/3) after s -> -s). The value is (Gamma(1/3) Gamma(2/3))^3 = (2pi/sqrt3)^3. So T(1) = K_B (2pi/sqrt3)^3 = 2pi/sqrt3.
- *Behaviour at u -> 0.* Let G_B(s) := K_B Gamma(1-s)^2 Gamma(1/3+s) Gamma(2/3+s), so that the integrand is -G_B(s) u^s/s^3 and G_B(0) = sqrt3/(2pi). The triple pole at s = 0 gives T(u) = (1/2)(G_B u^s)''(0) + o(1), with (log G_B)'(0) = -3 log 3 and (log G_B)''(0) = 2 psi'(1)+psi'(1/3)+psi'(2/3) = 5pi^2/3. So T(u) = (sqrt3/2pi)(l^2/2 + 5pi^2/6) + o(1).
- Hence t_c = (sqrt3/2pi)(l^2/2 + 5pi^2/6 - 4pi^2/3) + o(1) = (sqrt3/2pi)(l^2/2 - pi^2/2) + o(1), the same as t_cf, since g, g2 -> 0 as u -> 0. Both have theta = -w. Thus t_cf = t_c, and dt_c/dxi(0) = 1.

*Other branch.* On log y = l - i pi the same algebra gives, by direct substitution,

    (sqrt3/2pi)(I_1^{[2]} - i pi I_0^{[1]} - pi^2) = t_cf - sqrt3 pi - i sqrt3 (l + g).

l, g are real on (-1/27,0), so at xi = 0 the real part is -sqrt3 pi = -5.4414... != 0, and the printed formula does not vanish there.

This replaces the reliance of draft Lemma 1.2(c) on BS's printed expansion. In the coordinates of eq:tau this branch has the conifold at tau = 1, with A_cf = i sqrt3 (tau-1) A, while BS's description places it at tau = 0 (see s.5, item 7).

**C Step 7. Corollaries [P mod Theorem R_P2, Corollary SV_P2, Theorem VI, as in s.1.3].**

*Corollary G.* GW^cf_g = Sum_Gamma |Aut|^{-1} Sum_c T(Gamma,c) at c1 = C1cf. Each T(Gamma,c) is a meromorphic germ at t_c = 0 (draft Cor 2.8). Theorem R_P2 removes all polar terms of graphs with an edge, colouring by colouring. Corollary SV_P2: the single vertex is analytic at colours 1, 2, and at colour 0 its polar part is (3 eps1 eps2 E^3)^{g-1} Z[g,0] t_c^{2-2g} = (-3)^{g-1} beta_{2g-2} t_c^{2-2g}, using Theorem VI. With t_c = t_cf, this is the claim.

*Taylor coefficients in K.* wc, t_c and F_cf lie in Q[[xi]], and y^{-1} is regular at xi = 0. So an element of K[y^{+-1},Delta_c^{-1}][F_cf] has a Laurent expansion in xi over K, and t_c is a coordinate over Q. GW^cf_g is such an element, and t_c^{2-2g} lies in Q((xi)). So the difference has a Laurent expansion in xi over K, equivalently in t_c, and by Theorem R_P2 and Corollary SV_P2 it has no negative powers. F_cf and t_c converge near xi = 0, t_c is a local coordinate, and y^{-1} is regular there. So the difference is a finite expression that is meromorphic on a disc about xi = 0, and with no negative powers it is holomorphic.

*Corollary G_{k,h}.* By B5, eps1 eps2 GW_{g,d} is a symmetric homogeneous polynomial of degree 2g. For d = 0 this also holds, since lambda_g^2 = 0 removes the only 1/E. So eps1 eps2 F_g^{T_F} = Sum_j E^{2g-2j}(eps1 eps2)^j f_j(y), with f_j in Q[[y]]. The denominators in eps of the coefficients of F_g^{T_F} as a polynomial in F are monomials in eps1, eps2, E (from vertices and rows), so every rational eps with eps1 eps2 E != 0 is admissible. Solving at finitely many such eps (Vandermonde) puts f_j in Q[y^{+-1},Delta_c^{-1}][F], and uniqueness gives Cf(F_g^{T_F}) = Sum E^{2g-2j}(eps1 eps2)^{j-1} Cf(f_j). The Laurent coefficients of Cf(f_j) in t_c lie in Q, and the monomials E^{2k}(eps1 eps2)^{h-1}, k+h = g, are linearly independent over Q. Comparing at each polar order gives the statement. The limits F^cf_{0,h}, F^cf_{h,0} are coefficient extractions, so E = 0 and eps1 = 0 need not be points of T_F. With BS's d > 0 convention the left side changes by the constant GW^{T_F}_{g,0} in K, which affects only the O(1) term.

Identification, not used in the proof: by B5, eps1 eps2 GW^T_{g,d} is a polynomial, so it may be evaluated at alpha = 0, eps1 = -eps2 = eps. There e_T(R^1) becomes the non-equivariant Euler class of rank equal to the virtual dimension, and Lambda^vee(eps) Lambda^vee(-eps) = (-1)^g eps^{2g}. Hence the Q^d coefficient of F_{0,g}, d >= 1, is the usual local P2 invariant N_{g,d}. The d = 0 coefficient is the T_F constant and differs from other torus choices; this affects only the O(1) term.

**C Step 8. Comparison with BS [P by hand; comparison only].**
- rmk:cfgap display at g = 2: beta_2 = (7 eps1^4 + 10 eps1^2 eps2^2 + 7 eps2^4)/(5760 eps1 eps2), which gives exactly the printed coefficient.
- BS 5347, 5327-5332 are reproduced from Def 1.3 at xi -> 0 with signs. v_0 = -xi(1 + 2xi/3 + ...)/E gives Delta_0^{1/2} = -(3 eps1 eps2 E^3)^{1/2}(1/xi - 2/3 + ...), R^cf_00 = -xi h_w/((3 eps1 eps2 E)^{1/2}(E+z)), R^cf_10 = -h_w/(3 eps1 eps2 E)^{1/2}, R^cf_20 = (E/(3 eps1 eps2))^{1/2} h_w, using b(0) = 1, c(0) = -E, A_cf(0) = -1 from Prop 2.10 S6.
- beta conventions and BS 5299: the coefficient of (eps1 eps2)^{g-1} in (-3)^{g-1} beta_{2g-2} is 3^{g-1} B_{2g}/(2g(2g-2)), because beta_{2g-2}(eps,-eps) = eps^{2g-2} B_{2g}/(2g(2g-2)) and (eps1 eps2)^{g-1} = (-eps^2)^{g-1} there. The coefficient of E^{2g}/(eps1 eps2) is (-3)^{g-1} B_{2g}(1/2)/(2g(2g-1)(2g-2)), via x/(2 sinh(x/2)) = Sum B_n(1/2) x^n/n!.
- lem:cfpole is implied. The result here is sharper: every graph with an edge has no pole at all, and in the vertex graph only order 2g - 2 occurs.

**C Step 9. Reconciliation [P mod PPZ Props 2.8, 2.9; not an input].** Fpot_g(q_0) = Int T_b(h_w.1), with b = z(1 - E/(E-z)), and h_w.1 = R_h T_h 1. By Prop 2.9 T_b R_h = R_h T_{h^{-1} b}, and by Prop 2.8 the translations add: h^{-1}(z) b(z) + z(1 - h^{-1}(z)) = z(1 - h_w(-z) E/(E-z)). This is exactly BS's Cont_cf(d); the edge is Cont_cf(e). So eq:betagiv is Theorem VI. This resolves draft s.7.6 (B-c). prop:cfgraphsum is not used.

**Round-5 critic items, as hypotheses.** S1(a): true. S1(d): the literal claim is partly wrong (coloured terms do lie in BS's abstract G_{T_F}), but its substance is right: there is no map to functions without t. Handled in C Step 4. N3: resolved, with A irrelevant for point-free terms. The stored warning about a mechanism like h_nw: not needed.

---

## 5. Errata for BS and the draft

### 5.1 BS (BS_v2_2026-11-02/main.tex). None of the following is applied to that file.

1. **Lines 1277, 2446 versus eq:gtH (3799) and eq:F2ref (4217).** GW_g (line 1277) and the correlators of eq:gwprim (line 2446) are printed with d > 0, while the graph sum includes d = 0. F^T_2 at y^0 is non-zero, for example 19189/100800 at alpha = (1,2,5), eps = (3,7). The constants of eq:F2ref at y = 0 are 1/640 in F_{1,1} and 0 in F_{0,2}, F_{2,0}, which match E^2/640. Read GW_g with d >= 0, or read the identities up to the constant.
2. **Line 3537.** f_j^* must be the constant covector at Q = 0, not "at t".
3. **eq:normNLR (4327) with the printed Asym.** It holds with a constant c, c^2 = -1. With (det Hess W)^{-1/2} in place of (det Hess(-W))^{-1/2}, c = 1.
4. **Lines 4054-4056.** Regularity on the diagonal is asserted in one clause from properness of the T_F-fixed locus and lem: GW contravariance (1284), whose localisation is unspecified. The needed statement is localisation over the prime of T_F, which is not written. It is proved here another way: H^0 = 0 on the proper M, or e_3 for d = 0 (B5).
5. **Lines 4196-4198 (ansatz), 4345-4347 ("by continuity"), 3983 ("limiting version" of a rational torus) and 4056 ("in the closure of the set of rational weights").** Not needed, and replaced by B2, B4, B6.
6. **Lines 5163-5167.** Cf(f)(tau) = f(-1/(3tau)) is false as written. Under tau -> -1/(3tau), u = -27y -> 1 - u (from eq:tau and the eta transformation), so X -> 1/u, not X, contradicting line 5172. GW_g(-1/(3tau)) is regular at the conifold (it is GW_g near y = 0), contradicting thm:cfasym. F also fails, because theta does not commute with the Fricke map. On the generator A, one has A(-1/(3 tau')) = -A_cf with tau' = tau - 1, on the branch of Lemma (t_cf). The lines fixing X and Y_i (5159, 5174) are commented out. The substitution on generators, identity on coefficients, is the reading forced by 5172 and 5340.
7. **Line 5145.** t_cf needs the branch log y = log(-y) + i pi on (-1/27,0). On that branch eq:tau places the conifold at tau = 1, with A_cf = i sqrt3 (tau-1) A; it is T-equivalent to BS's tau = 0.
8. **Attribution of prop:algindABC (2854).** The independence of A, B, C is in BFGW s.4.1 (Lem 4.1 plus Zagier Prop 20). Lem 4.2 is its corollary, and it reads: S, X, I11 are algebraically independent. eq:ABCmod is called well known in BS and is cited to other papers in BFGW.
9. **Lines 4071-4076.** The passage argues regularity at y = 0 of the coefficients of the polynomial in F, but the given regularity is of the total, and F = O(y). So BS's Q[X^{+-1},F] does not follow as written. No counterexample is known (s.0.3).
10. **Printed sign of W_i on T_F (about 3382-3384).** Confirmed again by referee C; see draft s.9 item 1. The other items of draft s.9 stand.
11. **Lines 927-929.** The expansions of Pi_gamma are normalisation statements. They hold for an integral cycle only up to a constant prefactor, which BS does not state; at the conifold the prefactor is fixed by line 5145.
12. **Lines 2633, 2658.** MR2486673 treats only the fibre-rotating C^* (T_F); it does not give eq:J=I at generic alpha. CG07 line 2633 also has the torus on fibres only. The commented pointer Theorems 3.5, 3.7 at 2659 does not match the paper (concave case: Theorem 3.4, Prop 4.1, Cor 4.2). The correct citation at T_max is MR3414388, Theorem 31 (s.1.2).
13. **Line 2446 (eq:gwprim).** The double bracket is printed with d > 0, without 1/m! and without Q^d. eq:qprod, eq:fundsol, eq:Jred and the QDE need Sum_{d>=0, m>=0} (Q^d/m!) over stable (d, n+m). With the printed range the QDE is false.

### 5.2 Draft (LP2_REGULARITY_PROOF_draft.md): changes this file required; all eight applied in the draft at round 6c

1. s.0.2, s.7.6, s.8.1, s.8.3: Proposition 13.2_P2 and Corollary G are no longer [O]/conditional. Their status is that of s.0.2 here.
2. Prop 2.11: the proof is to be replaced by Theorem FB_B (a)-(d). The step "if delta divides Q some phi_{k,m} blows up" is to be deleted, not repaired. s.8.3 items 3 and 4 are closed. Item 4's algebraic independence is needed only for Cf to be well defined on functions, not for the identity, and on leg-free terms only the transcendence of F is needed.
3. Lemma 1.2(c): replaced by Lemma (t_cf), which proves the constant.
4. Prop 2.6: step 5 no longer rests on (C2) as asserted by BS, because (Normalisation) is A3.3. (C1) reduces to eq:J=I plus the QDE (A2.1-A2.4). (C2) and (Symplecticity) are discharged.
5. Thm 2.7(e): step 2 has the gap filled in B4. Step 5 is replaced by the single-variable identity of B4.
6. s.0.3: (C5) thm:haegen, lem:lho, eq:gtH at rational tori and BS 4056 are no longer used by Proposition 13.2_P2 (Stages A, B replace them). (C4) is used only for statements involving A.
7. s.1.3: "the Cf convention" should read "the reading of BS adopted as the definition" (Definition (Cf)).
8. s.0.4 and s.8.2: the bookkeeping of what has been refereed should follow s.1.3, s.1.4 here.

The corresponding conifold notes: conifold Prop 13.2's proof uses the Teleman numbers Theorem (2), 2.18, 2.19, Prop 8.7, with Prop 8.4, (8.6) for comparison and Prop 7.7 only in a remark, and the steps "Givental Prop (c)" and flatness as in A2.8-A2.9. A1.5, A1.7, A1.8, A2.9 are rank-free. A2.8 is rank-free except for the non-degeneracy Delta_i^g != Delta_k^g for some g, which is checked per target (conifold: (F8)). M_0 = 0 also follows with no non-degeneracy from A1.5 at (g,n) = (0,3): the coefficient of a_j b_i c_i in the degree-1 part of (R.omega)_{0,4}(a,b,c,H) is -Delta_i^{1/2}(lambda_i - lambda_j) r_ij, and the same coefficient in theta C_abc is Delta_i^{1/2} (W_Psi)_ji (the computation of the Sign check in A2.8, re-derived by hand in round 6d). The field-of-definition paragraph of the conifold Prop 13.2 was repaired in round 6d, after the major item of refCon in round6_scripts/round6c_result.json (passes on the repaired paragraph: hodge_vertex/lp2_explore/ROUND6_STATUS_2026-11-02.md). No statement of this file rests on that paragraph.

---

## Appendix. Evidence, stated as evidence only

Nothing in this appendix is a proof step. The arithmetic is exact (Fractions, Q(omega)) unless stated otherwise. Scripts are under hodge_vertex/lp2_explore/round6_scripts/<key>/, where <key> is FB_A, FB_B, FB_C (authors) or refFB_A, refFB_A2, refFB_B, refFB_C, refFB_C2 (referees, with code independent of the authors').

### A.1 Stage A: generic alpha against localisation

- **FB_A/g2d1.py.** Genus-2 graph sum (all 7 stable graphs), rows from eq:canR, at (alpha;eps) = (1,2,5;3,7), (2,-3,7;1,4), (1,3,8;5/2,-3/7), against localisation at d = 0, 1. Exact match in all three samples, for example 19189/100800 and 34231/1680 at the first. Poles in z cancel through y^2. Symplectic defects 0 for z^<=5, y^<=2. Control c_j -> -c_j fails at both degrees in all samples. Re-run by referee A with the same result. Near alpha = E/3, loc_d1 tends to the refined BPS d = 1 value: (0,1) gives 48241/1680 at (3,7) and 81/80 at (1,1).
- **FB_A/g2d2_limit.py.** Near the diagonal (delta = 1e-2, 1e-3, 1e-4, floating), F[2] + 6F[1] tends to the refined d = 2 value: -1596.42857 at (3,7), -54 at (1,1), 53.16494 at (5/2,-3/7). Not re-run by a referee.
- **FB_A/sp3d.py.** 3-dimensional formal stationary phase through z^3, y^2, 2 weight sets, 3 colours: omega(y,z) = (rho/rho_0)(y,z) prod_3 Theta(-z/a) exactly in all 6 cases. Control with Theta(+z/a) fails in all 6.
- **refFB_A/d0.py.** d = 0 value W_i^0 (p_3 - p_1^3)/17280 reproduces 19189/100800, 158713/290304, -10457963/29030400.
- **refFB_A2 (engine.py, rows.py, gsum.py, t1-t5.py).** Rows derived without BS's D, E: C^c_{11} solved from z -> infinity of S(H)|_{P_j}, through y^3; A^2 D = X and (D - X)/2 through y^3, one weight set. Graph sum against independent GP localisation, exact: g = 2, d <= 3 at (1,2,5;3,7), (2,-3,8;1,4); g = 3, d <= 2 at (1,2,5;3,7); g = 4, d <= 1 at (1,2,5;3,7) (this reaches z^7, that is a_4). With legs: <H^a>_1, <H^a>_2 through d = 2 at three weight sets. Hodge-vertex form equals point form at g = 2, 3 through y^2. Feynman Phi(alpha,Y;z) equals (Delta^{1/2} R~_{0i}) h_{(3 alpha_i, tau, tau')} through z^3, y^2, three colours, two weight sets; control h^{-1} fails. Constants: GW^{T_F}_{2,0} = E^2/640 and GW^{T_F}_{3,0} = E^4/161280 by exact interpolation at eps = (3,7), (1,1), (5/2,-3/7).
- **refFB_C2/gen.py.** Rows from the I-function at P_j, eq:RHb and (Normalisation), exact y-series. z-poles cancel; 2E(y) = D - X. Graph sum minus localisation is 0 at y^0..y^2 (g = 2, two alpha) and y^0..y^3 (g = 3). Controls (sign flip of c_j, dropping the -eps1, -eps2 weights, no h) fail, already at y^1. At g = 3, y^0 the sign flip is invisible, because dim Mbar_3 is even.
- **round6d/classlevel (own code, no project import; round6_scripts/round6d/classlevel/).** Defect M of A2.8 computed directly from the rows of A2.4-A2.6, with C^c_{1a} solved from the z = infinity expansion of the I-function. M = D M_rat D^{-1} with D = diag((W_j^0)^{1/2}) constant. M_rat = 0 over Q at 6 rational samples through z^5, y^3, and through z^7, y^5 at (1,2,5;3,7). M = 0 in F_p with the square roots as elements of F_p at 5 samples, p about 2^89, and through z^9, y^6 at (1,2,5;3,7), p about 2^127. Controls R^*, R(-z), H* -> -H*, theta -> -theta fail, R^* first at psi^1, y^1. h -> h^{-1} and no h give M = 0, as A2.9 says. Class level: GP localisation with psi classes pulled back from Mbar_{g,n} against the R.omega graph sum with legs, exact over Q: 150 comparisons per weight set at (g,n) = (1,1), (1,2), (2,1), (2,2) through d = 3 at (1,2,5;3,7), (1,3,8;5/2,-3/7), (3,7,-11/2;5/3,-2/7); n = 3 at g = 1, 2 through d = 2 at the first set; g = 3, n = 1 through d = 2 at the first set and d = 1 at the second; no mismatch. Also through d = 2 at two weight sets, no mismatch: lambda_1 insertions; Hodge-vertex form with legs; descendants at k <= 6 from ancestors through S_z of genus-0 GP. Divisor equation with the bubble term on both sides, n = 1 -> 2 through d = 3, n = 2 -> 3 through d = 2, two weight sets, no failure; the same equation without the bubble term fails. Controls that fail: legs at +psi, R^*, R(-z), h^{-1}, no h.
- **round6d/cite (round6_scripts/round6d/cite/scripts/lp_rows_check.py).** Comparison with Lho-Pandharipande, s.1.2: R_1 = (1-L^2)/18, R_2 = (1-24L-2L^2+25L^4)/648 compared exactly with eq:Ifun at H = lambda_0 = 1 through y^8; a control fails.
- **round6e/citeO1 (round6_scripts/round6e/citeO1/scripts/jI_check.py, output out_d4.txt).** Round 6e, exact fractions, finite order, evidence only. Direct graph sum of virtual localisation on Mbar_{0,n}(Z,d), n = 1, 2, d <= 4, independent of the I-function, at three weight points: alpha = (2/7,5/11,-3/13), (1/3,-4/5,11/2), and (1/7,3/5,11/13) with e_1(alpha) = eps1 + eps2. Results at all three: eq:J=I at P_0, P_1, P_2 through y^4 at a rational z; eq:Ian 1/z coefficients; eq:etaflat; S(1) = J/z and S(H) = theta_Q J from two-point sums; the QDE z theta_Q S(H^a) = Sum C^c_{1a} S(H^c), a = 0,1,2, through Q^4, with C from three-point invariants. Control: dropping the mirror map fails. GW_{0,d} eps1 eps2 = 3, -45/8, 244/9, -12333/64 at every point, and <H^2>_{0,1,d} eps1 eps2 / e_1(alpha) = -5, 87/4, -1544/9, 27671/16.

Sampling warning: alpha = (2,-3,7) has w_F + w_F' = 0 on fixed graphs with d >= 2 (alpha_1 - alpha_0 = -(alpha_2 - alpha_0)), so term-by-term localisation is singular there. It was used only for d <= 1. Avoid this sample at d >= 2.

### A.2 Stage B: rows, traces, T_F

- **FB_B/run1.py.** Phi(alpha*,Y;z) by 3-dimensional Wick expansion equals Theta(-z/E) S(v,z) through z^3 at (eps;y) = (3,7;1/5), (1,1;-2/7), (5/2,-3/7;3/11); control Theta(+z/E) fails. The genus-2 total as Tr (x) Tr in 𝒞 (x) 𝒞 at alpha*, with rows only from Phi and no y-series, equals eq:F2ref exactly at the same points and 𝖥 in {0, 1, -2/3, 5/3, 11/7}. Control Theta(-z/eps) fails at the one value tested per sample. (Pair) equals the pairing from explicit rows at A = 7/4, M, b <= 3.
- **FB_B/run2.py, re-run by referee B.** (a) At fixed y, 𝖥, total at alpha* + delta(1,2,-3) minus T_F value: 1.02e-3, 1.03e-5, 1.03e-7, 1.03e-9 for delta = 1e-1..1e-4: no blow-up at fixed y. (b) Generic alpha, y0 = 1e-25: Tr (x) Tr total minus GW_{2,0}, over y0, equals GW_{2,1} up to O(y0); control wrong h gives 0.3-0.45. (c) Along alpha = alpha* + s(1,2/3,-5/7), s = 1..13: 5th differences of eps1 eps2 GW_{2,1} and 8th differences of eps1 eps2 e_3 GW_{2,0} vanish; control without e_3 fails. T_F values at s = 0 by interpolation: GW_{2,1} = 48241/1680, 81/80, -13093109/26342400; GW_{2,0} = 5/32, 1/160, 841/125440. (d) y-series of eq:F2ref: y^0 = E^2/640, y^1 = refined d = 1, y^2 + 6y^1 = refined d = 2. (e) N(W)|_{y=0} = -1440 at alpha = (1,2,5); N(W) = E^9 y^2/Delta_c^2 at y = 1/5, -2/7, E = 10.
- **refFB_B/wick_ref.py (logs wick6, wick9).** Phi(alpha*,Y;z) = Theta(-z/E) S(v,z) exactly through z^9 at (eps;Y) = (3,7;-2/3), (2,-5/3;5/7). Phi(alpha,-alpha_i;z) = prod Theta(-z/a) through z^4 at three colours, alpha = (1,2,5).
- **refFB_B/phi_I.py, more_ref.py.** At alpha = (1,2,5), Phi(alpha,Y;z)/Phi(alpha,-alpha_0;z) against the stripped I-function, z <= 4 (z <= 5, three colours, K = 3), truncation K = 2..5 at y = 1.5e-8: differences fall with K as truncation only.
- **refFB_B/g2test.py, g3_ref.py, more_ref.py.** T_F total at g = 2 equals eq:F2ref exactly at (eps;t;𝖥) = (3,7;1/2;1), (1,1;-2/3;5/3), (5/2,-3/7;3/11;-2/3). Genus 3: 𝒯bar_3(F(y)) at t = 1e-3, 1e-4, 1e-5, 1e-7 agrees with GW_{3,0} + Q GW_{3,1} + Q^2 GW_{3,2} + Q^3 GW_{3,3} (refined BPS as external control), with stable residual (the next coefficient). GW_{3,0} = 125/2016 at (3,7), 707281/6195732480 at (5/2,-3/7), 1/10080 at (1,1). Flipping Theta(-z/E) or Theta(z/eps) fails. Hodge_v = h_v.1 exact at g = 2, 3, three- and five-element weight sets. GW^T_{g,0} = Sum_i (W_i^0)^{g-1} Int Hodge_{wt_i} at g = 2, 3: polynomiality along a line and interpolation equal the formula of (b) at (3,7), (5/2,-3/7).
- **refFB_C2/loc.py, cmp_loc.py, g4b.py.** Direct virtual localisation of K_P2 x C^2 (no mirror symmetry, no R-matrix) at generic alpha; exact polynomiality in alpha; interpolation to alpha = E/3. Tbar_g[F(y)] at t = +-1/50, +-1/500: (Tbar_g - Sum_{d<D} Q^d GW_{g,d})/Q^D tends from both signs of t to GW_{g,D}: g = 2, 3, D = 3 at (3,7), (1,2), and at (5/2,-3/7) for g = 2; for example, g = 3 at (3,7): 2827874.165, 2827876.721 at t = +-1/500 against 2827875.443. g = 4, D = 2 at (3,7), t = 1/500: -289820.888 against -289821.616. Controls (rows h_w^{-1}, rows without h_w, F + 1/9) fail.
- **refFB_C2/phi.py, phi_vs_I.py.** Own Wick expansion, exact through z^4. Theta(z/eps1) Theta(z/eps2) Phi at the root Y_j(y) equals the I-function row times h_{wt_j} up to O(y^6), scaling confirmed at y = 1e-5, 1e-7.

### A.3 Stage C: frame change, gap, t_cf

- **FB_C/hodge_lib.py (engine controls).** DVV: 1/1152, 1/82944, 7/240. Mumford recursion: lambda_1^3 = 1/2880 on Mbar_2, lambda_1^6 = 1/90720 on Mbar_3, 3! ch_3 = -1/5760 on Mbar_2. Graph counts 7, 42, and the stored 379 at g = 4.
- **FB_C/fc.py (log fc_g3.log).** Frame change at fixed t and numeric 𝖥, rows from the closed form, at (eps;t) = (3,7;1/5), (1,1;-2/7), (5/2,-3/7;3/11), 𝖥 in {0, 1, -2/3, 5/3, 11/7, -7/2, 2/9}. g = 2: 21/21 cases Tbar_2 = G_2 = stored F20 = eq:F2ref. g = 3: 21/21 cases Tbar_3 = G_3 = stored F30. (deg_𝖥 <= 6, so this is the polynomial identity at those points.) All graph values differ between frames (7 of 7, 42 of 42), so the identity is for the total only. Controls (Hodge vertices with rows h_w Phat; point vertices with rows h_w^{-1} Phat) fail in every case.
- **FB_C/fc4.py (log fc4.log).** g = 4, 379 graphs: point frame = Hodge frame = stored hodge_F40.m at (3,7;t = 1/5;𝖥 = 1) and (1,1;t = -2/7;𝖥 = -2/3). Value at the first: 73603363895647521165850828027/67984534801612800000000.
- **FB_C/gap_series.py.** Stored F20, F30, F40 at XX = F_cf as exact Laurent series in xi, against (-3)^{g-1} beta_{2g-2}/t_c^{2g-2}, at (3,7), (1,1), (5/2,-3/7), (1/3,1), (4/5,1): every polar coefficient of the difference is 0 at g = 2, 3, 4. Leading values reproduce gap_3 = -169/27216 at x0 = 1/3 and -1414951/224000000 at x0 = 4/5. Control XX -> XX + 1/9 fails in all cases; XX -> XX + xi/5 fails at g = 3, 4 and does not fail at g = 2 (s.A.4). First regular value, g = 2 at (3,7): -988747/4898880.
- **FB_C/cfgraph.py.** eq:betagiv graph sum = Z[g,0] = beta_{2g-2}/(-eps1 eps2 E^3)^{g-1} exactly at g = 2, 3, 4, (3,7), (1,1), (5/2,-3/7); with h(+z) the sum differs (g = 2, 3). Limit coefficients of beta at E = 0 and eps1 -> 0 equal BS 5299 at g = 2, 3, 4 (-1/80, 1/112, -3/160; -7/1920, -31/17920, -381/143360).
- **refFB_C2/tf.py.** Own graph enumeration. Tbar_2 = G_2 = eq:F2ref at random (eps,t,𝖥,c1) in Q(omega), c1-independence confirmed. Tbar_3 = G_3 at t = +-1/50, +-1/500 with F(y). Symplecticity with free c1, 𝖥 exact through z^8 on T_F. Controls fail.
- **refFB_C2/gapnum.py, gapnum4.py, g4.py (logs gap_g4.log, g4.log).** Exact evaluation of G_g[F_cf] at rational t = -1/3 + delta, no stored data. G_g - (-3)^{g-1} beta_{2g-2}/t_c^{2g-2} stays bounded and converges: g = 2 (delta = 1e-2..1e-5), g = 3 (1e-2..1e-4) at three eps pairs, and g = 4 at (3,7): 0.3303, 0.3858 at delta = 1e-2, 1e-3. Per colouring (g = 2, 3): max |T(Gamma,c)| over graphs with an edge is bounded; colour 1 vertex is bounded; colour 0 vertex minus polar term is bounded. Controls: F + 1/9 diverges; F + xi/5 diverges at g = 3; xi in place of t_c diverges.
- **FB_C/gap_series.py, refFB_C/tcf.py, refFB_C2/tcf.py (t_cf).** t_c minus BS's expression on log y = l + i pi: 2e-16, -7e-16, -4e-16 at y = -1/54, -1/81, -1/270 (double precision); below 7e-16 also at y = -1/40. 60-digit Decimal: agreement to 1e-58 at y = -1/54, -1/81, -1/40 (1e-38 at -1/270, truncation), and A_cf + wc = 1e-56. Other branch at y = -1/54: re = -4.6936, im = 6.6811, consistent with the identity of C Step 6.
- **Stored per-colouring data (RESULTS_2026-11-01.md, lp2_percolour_g2_out.txt, lp2_percolour_g3_out.txt).** See the draft's Appendix. Theorem R_P2 per colouring was tested at g = 2, 3 only.

- **round6d/cite (round6_scripts/round6d/cite/scripts/barnes_check.py, barnes_run2.py).** Barnes, parameters (2/3,2/3,2/3,4/3,5/3), trapezoid on Re sigma = -1/2, h = 1/120, 60 digits: LHS - (2 pi/sqrt 3)^3 = 8.1e-47; generic control (0.3,0.45,0.7,0.2): 8.6e-32; wrong e' = 2 fails (47.98 against 67.04). Stdlib decimal, independent of earlier checks [E].

### A.4 Limits and warnings

- **Genus 2 sees only F_cf(0) at the O(1) level.** The polar part of G_2 depends on F_cf only through its value at xi = 0, so the control F -> F + xi/5 stays bounded at g = 2. Genus 2 alone is not a test of F_cf beyond xi^0, and genus >= 3 is the first real test.
- **Rank-one sign controls do not discriminate at odd genus.** Int_{Mbar_g} Hodge_v is even in v when 3g-3 is even. So the control Theta(+z/v) versus Theta(-z/v) returned True at g = 3. Genus-3 sign tests must use the T_F total, where the flips do fail (refFB_B/g3_ref.py). Likewise, at g = 3, y^0 the sign flip of c_j is invisible.
- **R versus R^*.** R and R^* have the same R_1. Tests that see only first order in z cannot separate them. The genus >= 2 matches above do see z^2 and higher, but they remain evidence. LP/round6_scripts/round6c/refA/rstar.py, rstar2.py, rstar3.py (LP as in the draft's Appendix): exact divisor identity A Int_{Mbar_{g,n+1}} Omega(H^{n+1}) = theta_y Int_{Mbar_{g,n}} Omega(H^n) at alpha = (1,2,5), eps = (3,7). For R the defect is 0 at (g,n) = (0,4), (1,1), (2,0) through y^2 and at (3,0) through y^1. For R^* the defect is non-zero already at y^1: 1249/2268, 76487/17280, -28742463168731/35831808000, -406623543620925998705062409/62912004762894336000000. Evidence only. R^* differs from R only in class degree >= 2 of Omega. In round 6d, R^* passes every test at (g,n) = (1,1) and every test with total psi power k >= dim Mbar_{g,n} - 1 at (1,2), (2,1), (2,2). Those tests cannot separate R from R^*. Below that power it fails at (1,2) k = 0 (3 of 6), (2,1) k = 0, 1, 2 (2/3, 3/3, 3/3), (2,2) k = 0, 1, 2, 3 (3/6, 8/9, 15/15, 18/18), first at y^1.
- **Scale law.** Row and edge factors at different weight pairs differ only by eps1 eps2 and the scale E (draft Appendix caveat). Row and edge checks at n weight pairs are one test of the scale-free object plus a test of the scaling. Hodge-side quantities do depend on tau.
- **Coverage.** Nothing was checked to all orders by computation. Genus: g <= 3 for most items; g = 4 with one or two eps pairs and t points; g = 5 per colouring for graphs with a vertex of genus >= 2 at three eps pairs (round 6d); g = 6 for about 1000 of 5226 such graphs at one eps pair, with an unfinished log (round 6d); single-vertex polar part through g = 11 at (3,7) (round 6d). Nothing beyond these. Degree: d <= 3 (d <= 2 at g = 4; d <= 1 at g = 4 for generic alpha). z-order: z <= 3 (authors) and z <= 9 (referee B) for Wick rows on T_F; z <= 5 for rows at generic alpha against the I-function; z^7 for the constants a_k at generic alpha through g = 4. y-order: y^3 at generic alpha. Samples are small sets of rational weights. The refined BPS numbers are an external control, not used in any proof. Numeric regularity near T_F at fixed non-zero y was tested only at g = 2 (FB_B/run2.py (a)). Before round 6d the class-level identity Omega = R.omega was not tested beyond integrated numbers. Round 6d extends the rows at generic alpha to z^9, y^6 (one sample, mod p), and adds the class-level comparison with psi classes of s.A.1. It is still finite-order evidence.
- **Sensitivity to F_cf (round 6d, fresh code; round6_scripts/round6d/backward/, F_p with p = 2^127-1).** Replacing F_cf by F_cf + (1/9) xi^k in the edge graphs changes the sum at valuation k - (2g-3) at g = 2, 3, for k = 0..5. So the polar part sees F_cf through xi^{2g-4} at these genera, namely xi^0 at g = 2 and xi^2 at g = 3. At g = 5, with 356 edge graphs and 9484 coloured terms, xi^6 gives poles (3 colourings) and xi^7 gives none. Mode A gives poles in 2926 of 9484 terms. Evidence only.

### A.5 Disposition of the referee items

Every minor item of the five referee reports, by location of its fix. No report contained a major item.

| Report, item | Topic | Where |
|---|---|---|
| refFB_A 1 | colour isolation at n = 0 | A1.7, A3.1, A3.3 |
| refFB_A 2, 3 | ring and compatibility of the series | A1.3 (replaced by B5; defects recorded) |
| refFB_A 4 | log y; order of definitions of u_j | A2.4, A2.5, statement (A2) |
| refFB_A 5 | field for the kernel argument | A3.5 |
| refFB_A 6 | R versus R^*; PPZ (12); footnote 10 | A2.8 |
| refFB_A2 1 | log y; order of definitions | A2.4, A2.5 |
| refFB_A2 2 | "one power of mu" | A1.3 (moot after replacement) |
| refFB_A2 3 | footnote 10; sign from Def 7.1 | A2.8. Its clause "the first-order check fixes everything" is declined: see refFB_A 6 |
| refFB_A2 4 | fixed field of Aut | A1.8 |
| refFB_A2 5 | E_t := sigma(R) | A2.8 |
| refFB_A2 6 | rationality in z | A2.1 |
| refFB_A2 7 | dependency of A3.5 | A3.5, s.0.2 |
| refFB_A2 8 | values of the constants | Theorem FB_B(b), s.A.1 |
| refFB_A2 9 | sample (2,-3,7) | s.A.1 |
| refFB_B 1 | scope of the verdict | s.1.4 |
| refFB_B 2 | F transcendental, no modularity | s.3, proof of (e) |
| refFB_B 3 | transfer of Lemma 2.4 | A3.5, B1, B2 |
| refFB_B 4 | s-chart | B4 |
| refFB_B 5 | GP for T | B5 |
| refFB_B 6 | 𝖥 free; telescoping | B3, Theorem FB_B(a), C Step 0 |
| refFB_B 7 | integrally closed | s.3, proof of (e) |
| refFB_B 8 | wording of (d) | Theorem FB_B(d) |
| refFB_B 9 | rank-one control | s.A.4 |
| refFB_C 1 | BFGW Lem 4.2; attribution | s.1.2, C Step 4(i), s.5.1 item 8 |
| refFB_C 2 | elementary transcendence | s.3, proof of (e); C Step 4(i) |
| refFB_C 3 | Barnes; contour | C Step 6, s.1.2 |
| refFB_C 4 | erratum (a), last clause | s.5.1 item 6 |
| refFB_C 5 | reading of BS as definition | Definition (Cf) |
| refFB_C 6 | coloured Cf on lifts | Definition (Cf), C Step 4(ii) |
| refFB_C 7 | target of ev_A | Definition (Cf), C Step 4(i) |
| refFB_C 8 | tags; deletion; constant in eq:gtcf | C Steps 0, 2, 4(iii), 5 |
| refFB_C 9 | Mumford with markings | A3.2 |
| refFB_C2 1 | Barnes printed form; continuity | C Step 6, s.1.2 |
| refFB_C2 2 | transcendence through the logarithm | s.3, proof of (e) (the value of b) |
| refFB_C2 3 | sqrt(3 eps1 eps2) | Definition (Cf), C Step 4(iii) |
| refFB_C2 4 | identity on the other branch | C Step 6, Lemma (t_cf) |
| refFB_C2 5 | Taylor coefficients; admissible eps | C Step 7 |
| refFB_C2 6 | genus 2 limit of the evidence | s.A.4 |

Hand re-derivations made for this assembly, in addition to the referees': the colour isolation of A1.7, A3.3; the corrected u_j in A2.4; the Euler derivation of PPZ (12), including the different equation for R^*, and d_v Psi = -mu Psi from the charges; (K+D)^{2l} in A3.2; the transfer of A3.5; the s-chart power count and the SP5 factors of B4; the 𝖥-telescoping and (Pair) in B3; the monodromy argument in (e) and the valuation step; the Riccati equation; every formula of C Step 6, including both branches and the Barnes parameters in both forms; A_cf = i sqrt3 (tau-1) A from eq:tau; the Vandermonde and monomial step of C Step 7; the norm of W in B6.
