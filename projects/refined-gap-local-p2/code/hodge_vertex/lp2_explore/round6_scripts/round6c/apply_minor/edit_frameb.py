import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import Doc
d = Doc(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "FRAMEB_PROOF.md"))
r = d.rep

# Provenance, JSON keys
r("""**This assembled text has not had a referee pass of its own.**""",
  """**This assembled text had no referee pass of its own when it was written; the record of the round-6c passes is the marker ROUND6C-REFEREES in s.0.4 of the draft.** The JSON keys A, A2, B, C, C2 are refFB_A, refFB_A2, refFB_B, refFB_C, refFB_C2. Minor items of the round-6c referee reports are applied in place where they could be re-derived; A2.8 is the one place where a proof route changed.""")

# Rulings
r("""No finite-order computation is a proof step. Scope is the fibrewise torus only, per colouring.

Contents.""",
  """No finite-order computation is a proof step. Scope is the fibrewise torus only, per colouring. Stage A uses eq:J=I and the genus-0 QDE at generic alpha on T_max, not on T_F. Ruling Q1 names statements on T_F. eq:J=I at T_max is a genus-0 mirror theorem, with no TR, cited by BS to Coates, Commun. Math. Phys. 287 (2009) 1071-1108 (MR2486673), not opened here. This use is outside the literal wording of Q1 and is recorded here pending the user's confirmation.

Contents.""")

# Graph sums notation
r("""G_g[𝖥] and Tbar_g[𝖥] are defined in Theorem FB_B(d) and in s.4.""",
  """𝒯bar_g(F) is defined in Theorem FB_B(d). Tbar_g[𝖥] and G_g[𝖥] are defined in the Notation of s.4, and Tbar_g[𝖥] is the reduction of the polynomial 𝒯_g of Theorem FB_B(c).""")

# FB_A status
r("""[P mod BS eq:J=I and the genus-0 QDE; Teleman Theorem (2), Prop 7.7, Prop 8.7(ii); PPZ Defs 2.2, 2.13, Prop 2.14; Graber-Pandharipande; Mumford's theorem as in Faber-Pandharipande s.1.1; Faber-Pandharipande Theorem 4. Part A3.5 additionally mod draft SP1-SP5 and Lemmas 2.1, 2.3-2.5. Two referee passes on the report; the proof of Teleman's Prop 7.7 was read and not re-derived.]""",
  """[P mod BS eq:J=I and the genus-0 QDE, with eq:SHb; Teleman Theorem (2), Prop 8.7(ii), and 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss); PPZ Defs 2.2, 2.13, Prop 2.14; Graber-Pandharipande; the fundamental-class and splitting axioms of the virtual class, with Coates-Givental App. 1 (iv)-(vi) for the twist (Behrend-Manin not opened); Mumford's theorem as in Faber-Pandharipande s.1.1; Faber-Pandharipande Theorem 4. Part A3.5 additionally mod draft SP1-SP5 and Lemmas 2.1, 2.3-2.5. Two referee passes on the report form. A2.8 was rewritten in round 6c with a direct derivation in place of the report's route through Teleman's Prop 7.7; that derivation has had one hand re-derivation and no further pass.]""")

# A3 last sentence
r("""so R^Z = R^{K_P2} h_{(-eps1,-eps2)}, and Asym is R^{K_P2}_{0i}.""",
  """so R^Z = R^{K_P2} h_{(-eps1,-eps2)}, and Asym = c R^{K_P2}_{0i}, with R^{K_P2} taken in its own normalised frame, so that c = 1 in the (det Hess W)^{-1/2} convention. (R^Z_{0i} = (eps1 eps2)^{-1/2} R^{K_P2}_{0i} h_{(-eps1,-eps2)} and Theta(-z/eps) Theta(z/eps) = 1.)""")

# Graph sum wording
r("""(Graph sum, g >= 2.) In Kw[[Q]], colouring by colouring, with no square roots and no cube roots at generic alpha,""",
  """(Graph sum, g >= 2.) In Kw[[Q]] (each coloured term lies in Kw[[Q]], with no square roots and no cube roots at generic alpha; the identity is for the total),""")

# FB_B header
r("""Theorem A, Cor 2.2, Prop 2.10, Prop 2.12; Graber-Pandharipande. One referee pass, conditional on Stage A.]""",
  """Theorem A, Cor 2.2, Prop 2.10, Prop 2.12; Graber-Pandharipande; c(𝔼)c(𝔼^vee) = 1 (FP s.1.1), for the constant term. One referee pass, conditional on Stage A.]""")

# operator rename
r("""Qbar := Tbar Pbar_0 with Tbar := z theta + Y_i - (z/2) theta log Delta_i, Pbar_1 := Qbar/A, Pbar_2 := Delta_c A [Tbar Qbar - z F Qbar""",
  """Qbar := 𝐓bar Pbar_0 with 𝐓bar := z theta + Y_i - (z/2) theta log Delta_i, Pbar_1 := Qbar/A, Pbar_2 := Delta_c A [𝐓bar Qbar - z F Qbar""")
r("""pi(𝐓) = Tbar, since""", """pi(𝐓) = 𝐓bar, since""")
r("""After adjoining omega this equals Sum_c Tbar evaluated at t_v = omega^{c_v} t.""",
  """After adjoining omega this equals Sum_c of pi^{(x)n} T evaluated at t_v = omega^{c_v} t.""")
r("""[z^n] of Pbar_0, Qbar, Tbar Qbar are Laurent""", """[z^n] of Pbar_0, Qbar, 𝐓bar Qbar are Laurent""")

# (b)
r("""using (-3H)^{-1} = (H^2 + e_1 H + e_2)/(3 e_3) in H_T(P2)[e_3^{-1}]. Both lie""",
  """using (-3H)^{-1} = (H^2 + e_1 H + e_2)/(3 e_3) in H_T(P2)[e_3^{-1}], with H the equivariant class, H|P_i = -alpha_i. Both lie""")
r("""the T_F localisation definition. Constant term:

    GW^{T_F}_{g,0} = Int_{Mbar_g} [9 lambda_g lambda_{g-2} P(0) - 3 lambda_g lambda_{g-1} P'(0)],
    P(H) := Lambda^vee(E-3H) Lambda^vee(-eps1) Lambda^vee(-eps2)/((E-3H) eps1 eps2).

For g = 2 this is E^2/640 [P, by hand]. For g = 3, E^4/161280 at three weight pairs [E, s.A.2].""",
  """the T_F localisation definition. Constant term, where H' is the non-equivariant hyperplane class of P2 and Int_{P2} H'^2 = 1:

    GW^{T_F}_{g,0} = Int_{Mbar_g} [9 lambda_g lambda_{g-2} P(0) - 3 lambda_g lambda_{g-1} P'(0)],
    P(H') := Lambda^vee(E-3H') Lambda^vee(-eps1) Lambda^vee(-eps2)/((E-3H') eps1 eps2).

For g = 2 this is E^2/640 [P, by hand]. For g = 3 this is 9 E^4 Int lambda_1 lambda_2 lambda_3 = E^4/161280 [P, by hand, with FP Theorem 4; written in B5]. No closed form is claimed for g >= 4.""")

# (e)
r("""so the polynomial in 𝖥 is unique, and Cf(GW_g) := 𝒯bar_g(F_cf) is well defined.""",
  """so the polynomial in 𝖥 is unique, and Cf(F_g^{T_F}) := 𝒯bar_g(F_cf) is well defined.""")

# C1, C2 status
r("""[P mod Stages A, B and draft Prop 2.10, Mumford (through A3.2), PPZ Prop 2.14.] In L_omega[𝖥]""",
  """[P mod Stages A, B and draft Prop 2.10, Mumford (through A3.2), PPZ Defs 2.2, 2.5, 2.13, Prop 2.14.] In L_omega[𝖥]""")
r("""**Theorem FB_C (C2): dictionary (notation of s.4).** [P.] T'(z)""",
  """**Theorem FB_C (C2): dictionary (notation of s.4).** [P mod PPZ Defs 2.2, 2.5, 2.13 and draft Prop 2.10; at c1 = A also draft Prop 2.12.] T'(z)""")

# Prop 13.2_P2 status
r("""**Proposition 13.2_P2 (d >= 0 reading).** [P mod Stages A, B and draft Props 2.10, 2.12.]""",
  """**Proposition 13.2_P2 (d >= 0 reading).** [P mod Stages A, B; draft Props 2.10, 2.12; PPZ Defs 2.2, 2.5, 2.13, Prop 2.14 (C Steps 1, 2); and, for the clause F_cf = XX only, the logarithmic form in Lemma (t_cf) (DLMF 15.6.7, Gauss), not Barnes's second lemma.]""")

# Corollary G
r("""The latter three are P mod (MV) = LLZ eq. (1) and textbook facts, s.1.3.] For g >= 2,

    GW^cf_g = (-3)^{g-1}""",
  """The latter three are P mod (MV) = LLZ eq. (1) and textbook facts, s.1.3.] For g >= 2, with GW^cf_g := Cf(F_g^{T_F}), a germ at xi = 0, and with t_cf the analytic continuation to a neighbourhood of xi = 0 of BS's expression on (-1/27,0) with log y = log(-y) + i pi, which equals t_c by Lemma (t_cf),

    GW^cf_g = (-3)^{g-1}""")

# G_{k,h}
r("""This is (Conifold asymptotics) of BS def:dirint. In particular F^cf_{0,h}""",
  """This is (Conifold asymptotics) of BS def:dirint. Here Pi_{gamma_cf} is read through its displayed expansion: BS's Pi_{gamma_cf} solves L2b o theta = 0 and vanishes linearly at xi = 0 with leading coefficient 1 (BS 929, 937, 955-957). Since L2b has exponents (0,0) there, the solutions vanishing at xi = 0 with no logarithm are the multiples of t_c, so Pi_{gamma_cf} = t_c = t_cf. The coefficient 𝔟_{k,h} depends on this normalisation; the integrality wording of BS 925-927 is not used. F -> F_cf in def:dirint is ev_cf on the polynomial of Theorem FB_B(e), whose coefficients lie in Q[X^{+-1},(1-X)^{-1}] rather than BS's Q[X^{+-1}] (s.0.3). In particular F^cf_{0,h}""")

# s.0.3
r("""- the proof of (MV), Teleman's Prop 7.7 and eq:J=I, which stay cited.""",
  """- the proofs of (MV), of Teleman's Theorem (2), of Harer stability, Looijenga and Madsen-Weiss, and of eq:J=I, which stay cited.""")
r("""Dependency: Stage A -> Stage B -> Stage C (Proposition 13.2_P2), and Lemma (t_cf) is independent of A and B. Corollary G = Proposition 13.2_P2 + Lemma (t_cf) + Theorem R_P2 + Corollary SV_P2 + Theorem VI.""",
  """Dependency: Stage A -> Stage B -> Stage C (Proposition 13.2_P2). Lemma (t_cf) is independent of A and B. Its logarithmic form gives F_cf = XX in Proposition 13.2_P2; its constant (Barnes) gives t_c = t_cf in Corollary G. Corollary G = Proposition 13.2_P2 + Lemma (t_cf) + Corollary 2.8 + Theorem R_P2 + Corollary SV_P2 + Theorem VI.""")
r("""End-to-end inputs of Corollary G: eq:J=I and the genus-0 QDE; Teleman Theorem (2), Prop 7.7, Prop 8.7(ii); PPZ Defs 2.2, 2.5, 2.13, Props 2.14 (and 2.8, 2.9 only for the remark); Graber-Pandharipande; Mumford's""",
  """End-to-end inputs of Corollary G: eq:J=I and the genus-0 QDE, with eq:SHb; Teleman Theorem (2), Prop 8.7(ii), 2.18, 2.19; PPZ Defs 2.2, 2.5, 2.13, Prop 2.14 (and Props 2.8, 2.9 only for the remark); Graber-Pandharipande; the fundamental-class and splitting axioms of the virtual class, with Coates-Givental App. 1 (iv)-(vi) for the twist (Behrend-Manin not opened); Mumford's""")
r("""Every link has had at least one referee pass in its report form; s.1.4 lists""",
  """Every link has had at least one referee pass, in its report form or, for SP1-SP5 and Lemmas 2.1, 2.3-2.5, in an earlier round; s.1.4 lists""")

# s.1.1
r("""eq:J=I is cited by BS to MR2486673, which was **not opened**. It is not re-proved.""",
  """eq:J=I is cited by BS to MR2486673 (Coates, Commun. Math. Phys. 287 (2009) 1071-1108), which was **not opened**. It is not re-proved.""")
r("""BS's Asym (about 4245): definitions""", """BS's Asym (4256; Laplace expansion 4245-4252): definitions""")

# s.1.2 Teleman
r("""Proposition 8.7(i),(ii). **The proof of Proposition 7.7 was read, not re-derived; it is a sketch in places.** It is the one non-classical citation on which the graph sum rests. Journal numbering""",
  """Proposition 8.7(i),(ii); 2.17 (Harer stability), 2.18 (Looijenga), 2.19 (Madsen-Weiss), with the originals not opened. **The proof of Proposition 7.7 was read, not re-derived; it is a sketch in places. Since round 6c it is not an input:** A2.8 derives the flatness equation directly, and (7.7.a) appears there only as a comparison. Theorem (2) is a non-classical citation on which the graph sum rests; eq:J=I is a further citation. Journal numbering""")
r("""That the E = 0 and eps1 -> 0 extractions are those conjectures is BS's statement (line 5414), not checked here.""",
  """That the E = 0 and eps1 -> 0 extractions are those conjectures is BS's statement (line 5414), not checked here. The identification of F_{0,h} with the unrefined local P2 potential is in C Step 7; it is not used in the proof.""")

# s.1.3
r("""### 1.3 Draft inputs (current status, checked in the current draft s.0.2, s.7.6, s.8.1, s.8.2)

| Item""",
  """### 1.3 Draft inputs (current status, checked in the current draft s.0.2, s.7.6, s.8.1, s.8.2)

In s.1.3, s.1.4 "this round" means round 6b. Round 6a (D_rows, D_core, D_hodge, D_VI) re-derived Theorem A, Prop 2.10 S0-S6, Theorem LC_P2, s.7.5, Corollary SV_P2 and s.7.8 by hand; see draft s.0.4.

| Item""")
r("""Consequences tested through z^9 [E]. |""",
  """Consequences tested through z^9 [E]. Two referees in round 3, as F1-F5 of R3 Part C. Round-4 critic 1 re-derived Lemma 2.4 and did not re-derive the F1-F4 calculus in full. No pass in the present text at the time of writing. |""")
r("""| Thm 2.7(e), steps 1-3 | P mod (C1)-(C3) in the draft |""",
  """| Thm 2.7(e), steps 1-3 | P mod eq:J=I, genus-0 QDE as a whole theorem; steps 1-3 use only SP1-SP5, Lemma 2.4, Cor 2.2 and are not on the route except as FB:B4 |""")
r("""LLZ eq. (1) only, at mu = 1^a; LEMMA_A_PROOF.md Prop 2.1, Lemma 3.2, (4.2); textbook facts""",
  """LLZ eq. (1) only, at mu = 1^a; LEMMA_A_PROOF.md Prop 2.1, Lemma 3.2, (4.2) at k = 2; textbook facts""")

# s.1.4
r("""### 1.4 What the referee passes did not cover

- Proof of Teleman's Prop 7.7 and Theorem (2) (read only).""",
  """### 1.4 What the round-6b passes did not cover

- Proofs of Teleman's Theorem (2) and of 2.17-2.19 (Harer stability, Looijenga, Madsen-Weiss): used as stated, not re-derived. Prop 7.7 is not an input since round 6c; its proof was read, and its use is now a remark in A2.8.""")
r("""Stage A's formal-family use of Prop 7.7 was examined by referees A, A2 only.""",
  """Stage A's former formal-family use of Prop 7.7 was examined by referees A, A2 only. The direct derivation now in A2.8 has had no further pass. Referee C also took Theorems LC_P2, R_P2, SV_P2, VI as proved and read Stages A, B only at their interface.""")
r("""- The class-level identity Omega = R.omega beyond integrated numbers (n <= 1, no psi insertions at d > 0) was not tested numerically and rests on Teleman.""",
  """- The class-level identity Omega = R.omega beyond integrated numbers (n <= 1 in the round-6b runs, no psi insertions at d > 0) was not tested numerically and rests on Teleman. The R versus R^* test of s.A.4 reaches n = 5 at genus 0 and n = 2 at genus 1, still with no psi insertions.""")

# A1.5
r("""**A1.5 Divisor [P].**""", """**A1.5 Divisor [P mod the fundamental-class axiom].**""")

# A1.7
r("""So ch_{2k-1}(𝔼) != 0 on Mbar_{k,1} by A3.3 (for k = 1, Mbar_{1,1} is used directly). Hence r_i = 0.""",
  """So ch_{2k-1}(𝔼) != 0 on Mbar_{k,1} by the first paragraph of A3.3, which uses only c(𝔼)c(𝔼^vee) = 1 and FP Theorem 4, not A1.7, A1.8 or A3.1 (for k = 1, Mbar_{1,1} is used directly). Hence r_i = 0.""")

# A1.8
r("""so it embeds in ℂ. Teleman Theorem (2)""",
  """so it embeds in ℂ. Since H^*(Mbar_{g,n}) is finite dimensional, iota(Omega) is a CohFT with flat unit in H^*(Mbar_{g,n};Q) (x) ℂ, semisimple because Omega is semisimple over F_0. PPZ Defs 2.2, 2.5, 2.13 and Prop 2.14 are identities of graph sums, and their proofs hold over any field of characteristic 0. Teleman Theorem (2)""")

# A2.1 row 2
r("""This covers the third QDE row used in A2.4. So all these identities may be expanded at z = 0.""",
  """This covers the third QDE row used in A2.4. So all these identities may be expanded at z = 0.

*Row 2, with no lem:modrel.* With s = 1/z, eq:Ifun gives S(1) = I/z = 1 + s H I_0^{[1]} + s^2 (H I_1^{[1]} + H^2 I_1^{[2]}) + O(s^3) (draft Prop 2.12, Step 0). So S(H) = z theta_Q S(1) = H + s (e_1 (G/A) H + (P_I/A) H^2) + O(s^2), with G, P_I as in that Step 0. Since S(H^c) = H^c + O(1/z), comparing the z^0 terms of z theta_Q S(H) = C^0_11 S(1) + C^1_11 S(H) + C^2_11 S(H^2) gives C^0_11 = 0, C^1_11 = e_1 E(y)/A, C^2_11 = D/A, by eq:ADEFXdef alone. D(0) = 1, so S(H^2) = (z theta_y S(H) - e_1 E(y) S(H))/D. This is the content of draft Remark 2.13(a), used here as a proof step.""")

# A2.4 rename
r("""Put w := exp(-(u_j + alpha_j log y)/z) s^_{.j} = exp(-U_j/z) s^_{.j}, in K~w((z))[[Q]].""",
  """Put w^ := exp(-(u_j + alpha_j log y)/z) s^_{.j} = exp(-U_j/z) s^_{.j}, in K~w((z))[[Q]].""")
r("""Then w = e^{-u_j/z} s_{.j}, and from the QDE z theta_Q w = (C - lambda_j) w, with w|_{Q=0} = Psi_{.j}(0).""",
  """Then w^ = e^{-u_j/z} s_{.j}, and from the QDE z theta_Q w^ = (C - lambda_j) w^, with w^|_{Q=0} = Psi_{.j}(0).""")
r("""So w is unique and w = v. This proves""", """So w^ is unique and w^ = v. This proves""")

# A2.6
r("""Inverting the matrix identity and removing the exponentials (i = j cancels; i != j gives zero) yields Sum R~_{ai}(-z) eta^{ab} R~_{bj}(z) = delta_ij.""",
  """Write s = s^ diag(y^{-alpha_k/z}). Then G = s^(z) s^(-z)^T with s^ in GL_3(K~w(z)[[Q]]), since s^(Q=0) is invertible. Inverting G = eta gives s^(-z)^T eta^{-1} s^(z) = 1. Conjugating by exp(U_k/z), that is R~_{ai}(-z) = e^{U_i/z} s^_{ai}(-z) and R~_{bj}(z) = e^{-U_j/z} s^_{bj}(z), yields Sum R~_{ai}(-z) eta^{ab} R~_{bj}(z) = delta_ij.""")

# A2.8
i0 = d.s.index("""**A2.8 Flatness of Teleman's R [C, P].**""")
i1 = d.s.index("""*Sign check, first order only [P].*""")
new_a28 = """**A2.8 Flatness of R [P mod Teleman Theorem (2) through A1.8, and Teleman 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss)].**

Rewritten in round 6c. This derivation does not use Teleman's Prop 7.7. It has had one hand re-derivation and no further pass. In A2.8 only, D denotes the diagonal matrix below, not BS's D(y).

*Set-up.* Work over F_0 with the derivation theta := theta_Q, using Omega = R.omega (A1.8) and A1.5. Put B(psi) := R^{-1}(psi) and a^i(psi;x) := f_i^*(B(psi)x) for x in V over Kw, so theta x = 0. On M_{g,n} (smooth curves) only the edgeless graph survives:

    Omega_{g,n}(x_1,...,x_n) = Sum_i D_i N_i(kappa) prod_l a^i(psi_l;x_l),

with D_i := Delta_i^{(2g-2+n)/2} and N_i(kappa) := Sum_m (1/m!) p_{m*} prod_k (Delta_i^{1/2} f_i^* T(psi_{n+k})), a polynomial in the kappa_j with constant term 1.

*Pushforward.* Let p: C -> M_{g,n} be the restriction of Mbar_{g,n+1} -> Mbar_{g,n}, with sections delta_l. On C: psi_l = p^*psi_l + delta_l, psi_{n+1} delta_l = 0, delta_l delta_l' = 0 for l != l', kappa_j = p^*kappa_j + psi_{n+1}^j, and the only graphs with an edge are the bubbles (0,3) carrying legs l, n+1. So B(psi_l) = B(p^*psi_l) + delta_l (B(psi_l) - 1)/psi_l, and the bubble puts (1 - B(psi_l))/psi_l applied to H*x_l on leg l of the genus-g vertex (its edge factor at psi_h = 0 is (eta^{-1} - (1 (x) B(psi_{h'})) eta^{-1})/psi_{h'}). Since Delta_i^{(2g-1+n)/2} a^i b^i = Delta_i^{(2g-2+n)/2} f_i^*(a*b), the two delta_l terms add up to (1/psi_l)[H*,B(psi_l)] x_l on leg l. Hence

    p_* Omega_{g,n+1}(x,H) = Sum_i D_i N_i Sum_l f_i^*((1/psi_l)[H*,B(psi_l)]x_l) prod_{l' != l} a^i(psi_l';x_l') + Sum_i G_i(kappa) prod_l a^i(psi_l;x_l),

with G_i a polynomial in kappa.

*The defect M.* Define M(psi) in Mat_3(F_0[[psi]]) by

    theta a^i(psi;x) - f_i^*((1/psi)[H*,B(psi)]x) = Sum_k M(psi)^i_k a^k(psi;x),

which is possible because x -> a(psi;x) is invertible. A1.5 at n = 2, where D_i = Delta_i^g, reads

    Sum_i D_i N_i [(M(psi_1)a_1)^i a_2^i + a_1^i (M(psi_2)a_2)^i] = Sum_i Gtilde_i(kappa) a_1^i a_2^i,   Gtilde_i := G_i - theta(D_i N_i).

*Stability.* By Teleman 2.18, 2.19, monomials in kappa, psi_1, psi_2 are independent in H^*(M_{g,2};F_0) in total degree < (g-1)/2. Take the coefficient of kappa^0 psi_2^0, where N_i = 1 and a_2 = (f_i^* x_2)_i is arbitrary, and let x_1 be arbitrary. With D = diag(D_i):

    D M(psi) + M_0^T D = Gtilde(0), diagonal, modulo that degree.

Since g is arbitrary and M does not depend on g, M_n = 0 for n >= 1. Further (M_0)_ik = -(W_Psi)_ik + ([Lambda, R_1])_ik in the f-frame, and both terms are skew (R_1 is symmetric by symplecticity). So (Delta_i^g - Delta_k^g)(M_0)_ik = 0, and Delta_i^g != Delta_k^g because W_i^0/W_k^0 is a non-constant element of Kw, not a root of unity. Hence M = 0, that is

    theta(Psi^{-1} B(psi)) = Psi^{-1} (1/psi)[H*,B(psi)].

*The equation.* Taking eta-adjoints, with B(psi)^dagger = R(-psi) and Psi^dagger = Psi^{-1}, and putting z = -psi gives

    z theta_Q(R(z) Psi) = (H* R(z) - R(z) H*) Psi.

This is an equation for R with legs R^{-1}(psi), not for R^*. With H* Psi = Psi Lambda and theta_Q U = Lambda, it is equivalent to z theta_Q S = H* S for S = R Psi e^{U/z}, which is Givental (3).

*Remark: comparison with Teleman's Prop 7.7 [C; not a proof step].* Put E_t := sigma(R), over F_0[[t]]. Applying sigma to coefficients gives sigma(R).sigma(omega) = sigma(Omega). Apply the embedding iota: F_0 -> ℂ of A1.8 coefficientwise in t. Then iota(sigma(Omega)) is Def 7.1 for Z = iota(Omega) and u = -tH over ℂ[[t]], iota(sigma(R)) is its data, and the coefficient of t^0 in (7.7.a) is iota of the identity over F_0 stated above. theta_Q is never extended to ℂ. Only the direction H is used, because the integrand p_* Omega_{g,n+1}(x,H) belongs to the same family. Def 7.1 has the factor (-1)^m, and the divisor equation used is the true one, A1.5; the sign printed in Teleman Ex. 7.17 is not relied on. By A1.7 E_t is the unique such data. (7.7.a) along the line reads d_t(𝐄 Psi) Psi^{-1} = -[𝐄, H*_t]/z, and at t = 0 this is the equation above, which is also (7.7.a) along v = -H. The legs of (7.6) carry 𝐄_u^{-1}(psi_i), and the legs of PPZ Def 2.2 carry R^{-1}(psi_i), so in this comparison Teleman's 𝐄 is PPZ's R, with the same z; this is a reading of (7.6) against Def 2.2. Teleman's Prop 8.2(ii) is consistent with this. Theorem (2) is a non-classical citation from Teleman on which the graph sum rests; eq:J=I is a further citation.

"""
d.s = d.s[:i0] + new_a28 + d.s[i1:]
d.n += 1

# A3.1
r("""**A3.1 Reduction mod Q [P, GP].** All of R, Delta_i, e_i, eta^{-1} are regular at Q = 0, and the graph sum""",
  """**A3.1 Reduction mod Q [P, GP].** R (A2.9), R^{-1} (symplecticity), e_i, Delta_i^{+-1/2} and eta^{-1} are regular at Q = 0 (Delta_i(0) = W_i^0 != 0), and the graph sum""")
r("""By A2.3, R|_{Q=0} = diag(s_i), s_i := exp""",
  """By A2.4 (R~|_{Q=0} = Psi|_{Q=0}) and A2.9, R|_{Q=0} = diag(s_i) on the f_i(0), s_i := exp""")

# A3.3
r("""Taking t^{2g-2} and using lambda_g^2 = 0 gives (2g-1)!""", """Taking t^{2g-2} gives (2g-1)!""")

# Discharged
r("""**Discharged for the draft.** (C2) is now [P] (A3.3).""",
  """**Discharged for the draft.** (C2): (Normalisation) is now a theorem, A3.3, modulo the inputs of A2.8, A2.9, A3.2, A3.3.""")

# A3.5
r("""the theta-constants are the coefficient field, and W is a unit.""",
  """the theta-constants are the coefficient field, and Y, W are units, since Y_i(0) = -alpha_i and W_i(0) != 0 (the algebraic form of Lemma 2.5, stated after its proof in the draft).""")
r("""By Lemma 2.3, 𝒜 = c_0 W^{-1/2} Phi(alpha,Y;z), with [z^n] Phi in Q[alpha,Y,W^{-1}].""",
  """By Lemma 2.3, 𝒜 = c_0 W^{-1/2} Phi(alpha,Y;z), with c_0^2 = -1/3 and [z^n] Phi in Q[alpha,Y,W^{-1}].""")
r("""So the z-coefficients of W^{1/2} 𝒫 (W^{-1/2} Phi) lie in Q[alpha,Y,W^{-1}]. Lhat is linear in y with coefficient 27 Y^3, so y = -prod_j (Y+alpha_j)/(27 Y^3), and the fraction""",
  """So the z-coefficients of W^{1/2} 𝒫 (W^{-1/2} Phi) lie in Q[alpha,Y^{+-1},W^{-1}], because y = -prod_j (Y+alpha_j)/(27 Y^3) (Lhat is linear in y with coefficient 27 Y^3), and the fraction""")
r("""Applying iota_i (Y -> Y_i in Kw[[y]], with W(Y_i) a unit, commuting with theta) gives 𝒫_i (W(Y_i)^{-1/2} iota_i(Phi)) = 0 in K^[[y]][[z]].""",
  """Applying iota_i (Y -> Y_i in Kw[[y]], with Y_i and W(Y_i) units, since Y_i(0) = -alpha_i and W(-alpha_i) != 0, commuting with theta) gives 𝒫_i ((3W(Y_i))^{-1/2} iota_i(Phi)) = 0 in K^[[y]][[z]]. Here (3W(Y_i))^{1/2} lies in K^[[y]], because 3W(-alpha_i) is the product of the three a in the definition of K^. sqrt 3 itself is not adjoined.""")
r("""*Kernel argument.* Put A_i := c_0 W(Y_i)^{-1/2} iota_i(Phi).""",
  """*Kernel argument.* Put A_i := c (3W(Y_i))^{-1/2} iota_i(Phi), which equals c_0 W(Y_i)^{-1/2} iota_i(Phi) with c_0 = c/sqrt 3.""")
r("""*Value at y = 0.* Set q_i = y p. The jets of the phase at the critical point in log coordinates are those of W_0 = (w' - 3 alpha_i log w') + Sum_{l != i}(q_l - (alpha_l - alpha_i) log q_l), with w' = (p q_j q_k)^{1/3}, dOmega = dlog w' dlog q_j dlog q_k (SP1, SP4). Each factor""",
  """*Value at y = 0.* Evaluate Phi at Y = -alpha_i, an algebraic substitution with no limit, since W(-alpha_i) = alpha_i (alpha_j - alpha_i)(alpha_k - alpha_i) != 0. There H = diag(0,a,b) + (alpha_i/3)J and T_m = a xi_j^m/m! + b xi_k^m/m! + 3 alpha_i s^m/m!, with a = alpha_j - alpha_i, b = alpha_k - alpha_i, s = (xi_i + xi_j + xi_k)/3. In the linear coordinates (s, xi_j, xi_k) the phase splits, so SP1 (linear case) and SP4 factorise Phi into three one-dimensional expansions. Each factor""")

# B1
r("""So the ring is 𝒞 (a local etale O-algebra), no roots""",
  """So the ring is 𝒞 (a local etale O-algebra, because Lhat'(Y) = -W/Y is a unit), no roots""")

# B2
r("""- Stage A, in the form written in A3.5: c_0 W(Y_i)^{-1/2} iota_i(Phi) = c N R_{0i} in K^[[y]][[z]].""",
  """- Stage A, in the form written in A3.5: c (3W(Y_i))^{-1/2} iota_i(Phi) = c N R_{0i} in K^[[y]][[z]], where c (3W(Y_i))^{-1/2} = c_0 W(Y_i)^{-1/2} with c_0 = c/sqrt 3.""")

# B3
r("""a tensor in 𝒞 (x)_O 𝒞[𝖥], with no idempotent needed. For a loop,""",
  """a tensor in 𝒞 (x)_O 𝒞[𝖥], with no idempotent needed. Fix an orientation of every edge to define T_{Gamma,dec}. The map 𝒞 (x)_O 𝒞 -> Kw((y))^9 given by the iota_i (x) iota_j is injective, since 𝒞 (x)_O 𝒞 is free over O and O embeds in Kw((y)). Hence the tensor of E(a,b) equals the flip of the tensor of E(b,a) in 𝒞 (x)_O 𝒞[𝖥], T_{Gamma,dec} does not depend on the orientation, and this passes to L (x)_{K_y} L under pi. For a loop,""")

# B5 constant term
r("""*Constant term.* Mod H^3, c_top(𝔼^vee box T_P2) = -3H lambda_g lambda_{g-1} + 9 H^2 lambda_g lambda_{g-2}, using x1+x2 = 3H, x1 x2 = 3H^2, lambda_g^2 = 0, lambda_{g-1}^2 = 2 lambda_g lambda_{g-2} (from c(𝔼) c(𝔼^vee) = 1 [C]). In P(H) the term with 1/(E-3H) carries lambda_g, so it dies against lambda_g in the prefactor: eps1 eps2 GW^{T_F}_{g,0} is a polynomial, with no 1/E. For g = 2 this gives 9 E^2 Int lambda_1 lambda_2 = E^2/640.""",
  """*Constant term.* Here H' is the non-equivariant hyperplane class. Mod H'^3, c_top(𝔼^vee box T_P2) = -3H' lambda_g lambda_{g-1} + 9 H'^2 lambda_g lambda_{g-2}, using x1+x2 = 3H', x1 x2 = 3H'^2, lambda_g^2 = 0, lambda_{g-1}^2 = 2 lambda_g lambda_{g-2} (from c(𝔼) c(𝔼^vee) = 1 [C]). In P(H') the term with 1/(E-3H') carries lambda_g, so it dies against lambda_g in the prefactor: eps1 eps2 GW^{T_F}_{g,0} is a polynomial, with no 1/E. For g = 2 this gives 9 E^2 Int lambda_1 lambda_2 = E^2/640. For g = 3, put u = E - 3H', p = eps1 eps2. Mod lambda_3 and with lambda_1^2 = 2 lambda_2, P(H') = (u^2 - lambda_1 u + lambda_2)(p^2 + lambda_1 p E + lambda_2 E^2 + ...), so [P(0)]_2 = lambda_2 (E^2 - p)^2 and [P'(0)]_1 = -3 lambda_1 (2 E^2 p - p^2). Hence GW^{T_F}_{3,0} = 9 I [(E^2 - p)^2 + 2 E^2 p - p^2] = 9 I E^4, with I = Int lambda_1 lambda_2 lambda_3 = 1/1451520 (FP Theorem 4).""")

# B6
r("""- By B1 and B4 the reduction of 𝒯_g is the T_F colour sum. This is (d).""",
  """- By B1 and B4 the reduction of 𝒯_g is the T_F colour sum, also coefficientwise in the indeterminate 𝖥, so the reduction of the polynomial is Tbar_g[𝖥] of s.4. This is (d).""")

# C Step 2
r("""Tbar_g[𝖥] = G_g[𝖥] in L'. Both sides lie in L_omega[𝖥] by Remark (signs) of the draft.""",
  """Tbar_g[𝖥] = G_g[𝖥] in L'. Both sides lie in L_omega[𝖥]: by (Pair) of C Step 0 the edges do not contain c1, Dil uses only S, and by Remark (signs) of the draft kappa_0 does not survive.""")

# C Step 4(ii)
r("""ev_cf embeds L~ by the branch at xi = 0, which defines colour 0. For real t""",
  """ev_cf embeds L~ by the branch at xi = 0, which defines colour 0. This is the branch of draft s.0.1 and Lemma 1.2(b), with t_i = omega^i t and t -> -1/3. There v_0 = (1+3t)/(E t) = -xi/E + O(xi^2), v_1 and v_2 are non-zero at xi = 0, and Y_0 is the expression at BS 5340. So the colour called 0 under ev_cf is the singular colour of Theorems LC_P2, R_P2 and Corollary SV_P2. The identification is needed only to match per-colour statements, since the colour-summed total is independent of it. For real t""")

# C Step 5
r("""**C Step 5. Proposition 13.2_P2 [P mod Stages A, B].**""",
  """**C Step 5. Proposition 13.2_P2 [P mod the inputs in its label, s.0.2].**""")
r("""With the Lemma below, A_cf = C1cf and F_cf = XX, so""",
  """With the logarithmic form in the Lemma below (C Step 6, no Barnes), A_cf = C1cf and F_cf = XX, so""")

# C Step 7
r("""The difference in Corollary G is such an element, up to the constant in K, and it has no negative powers. So it is holomorphic with Taylor coefficients in K.""",
  """GW^cf_g is such an element, up to the constant, and t_c^{2-2g} lies in Q((xi)). So the difference has a Laurent expansion in xi over K, equivalently in t_c, and by Theorem R_P2 and Corollary SV_P2 it has no negative powers. F_cf and t_c converge near xi = 0, t_c is a local coordinate, and y^{-1} is regular there. So the difference is a finite expression that is meromorphic on a disc about xi = 0, and with no negative powers it is holomorphic.""")
r("""and uniqueness gives Cf(GW_g) = Sum E^{2g-2j}(eps1 eps2)^{j-1} Cf(f_j).""",
  """and uniqueness gives Cf(F_g^{T_F}) = Sum E^{2g-2j}(eps1 eps2)^{j-1} Cf(f_j).""")
r("""The limits F^cf_{0,h}, F^cf_{h,0} are coefficient extractions, so E = 0 and eps1 = 0 need not be points of T_F.

**C Step 8.""",
  """The limits F^cf_{0,h}, F^cf_{h,0} are coefficient extractions, so E = 0 and eps1 = 0 need not be points of T_F. With BS's d > 0 convention the left side changes by the constant GW^{T_F}_{g,0} in K, which affects only the O(1) term.

Identification, not used in the proof: by B5, eps1 eps2 GW^T_{g,d} is a polynomial, so it may be evaluated at alpha = 0, eps1 = -eps2 = eps. There e_T(R^1) becomes the non-equivariant Euler class of rank equal to the virtual dimension, and Lambda^vee(eps) Lambda^vee(-eps) = (-1)^g eps^{2g}. Hence the Q^d coefficient of F_{0,g}, d >= 1, is the usual local P2 invariant N_{g,d}. The d = 0 coefficient is the T_F constant and differs from other torus choices; this affects only the O(1) term.

**C Step 8.""")

# C Step 8
r("""- beta conventions and BS 5299: E = 0 gives 3^{g-1} B_{2g}/(2g(2g-2)). eps1 -> 0 gives (-3)^{g-1} B_{2g}(1/2)/(2g(2g-1)(2g-2)), via x/(2 sinh(x/2)) = Sum B_n(1/2) x^n/n!.""",
  """- beta conventions and BS 5299: the coefficient of (eps1 eps2)^{g-1} in (-3)^{g-1} beta_{2g-2} is 3^{g-1} B_{2g}/(2g(2g-2)), because beta_{2g-2}(eps,-eps) = eps^{2g-2} B_{2g}/(2g(2g-2)) and (eps1 eps2)^{g-1} = (-eps^2)^{g-1} there. The coefficient of E^{2g}/(eps1 eps2) is (-3)^{g-1} B_{2g}(1/2)/(2g(2g-1)(2g-2)), via x/(2 sinh(x/2)) = Sum B_n(1/2) x^n/n!.""")

# s.5.1 item 4, 5, 11
r("""4. **Lines 4055-4057.** The assertion that regularity on the diagonal follows "by properness" is incomplete: a proper fixed locus alone does not give regularity at weights where the T-fixed loci change. What makes it true is H^0 = 0 on the proper M, or e_3 for d = 0 (B5).""",
  """4. **Lines 4054-4056.** Regularity on the diagonal is asserted in one clause from properness of the T_F-fixed locus and lem: GW contravariance (1284), whose localisation is unspecified. The needed statement is localisation over the prime of T_F, which is not written. It is proved here another way: H^0 = 0 on the proper M, or e_3 for d = 0 (B5).""")
r("""4345-4347 ("by continuity"), 3983 ("closure of rational tori").** Not needed""",
  """4345-4347 ("by continuity"), 3983 ("limiting version" of a rational torus) and 4056 ("in the closure of the set of rational weights").** Not needed""")
r("""10. **Printed sign of W_i on T_F (about 3382-3384).** Confirmed again by referee C; see draft s.9 item 1. The other items of draft s.9 stand.""",
  """10. **Printed sign of W_i on T_F (about 3382-3384).** Confirmed again by referee C; see draft s.9 item 1. The other items of draft s.9 stand.
11. **Lines 925-927.** The expansions of Pi_gamma are normalisation statements. They hold for an integral cycle only up to a constant prefactor, which BS does not state; at the conifold the prefactor is fixed by line 5145.""")

# s.5.2
r("""### 5.2 Draft (LP2_REGULARITY_PROOF_draft.md): changes this file requires""",
  """### 5.2 Draft (LP2_REGULARITY_PROOF_draft.md): changes this file required; all eight applied in the draft at round 6c""")
r("""and the steps "Givental Prop (c)" and flatness as in A2.8-A2.9.""",
  """and the steps "Givental Prop (c)" and flatness as in A2.8-A2.9. A1.5, A1.7, A1.8, A2.8, A2.9 are rank-free. The field-of-definition paragraph of the conifold Prop 13.2 carries an OPEN mark from round 6c (sign of u); no statement of this file rests on that paragraph.""")

# A.4 R vs R^*
r("""The genus >= 2 matches above do see z^2 and higher, but they remain evidence.""",
  """The genus >= 2 matches above do see z^2 and higher, but they remain evidence. SP/lp2_round6/round6c/refA/rstar.py, rstar2.py, rstar3.py (SP as in the draft's Appendix; scratchpad): exact divisor identity A Int_{Mbar_{g,n+1}} Omega(H^{n+1}) = theta_y Int_{Mbar_{g,n}} Omega(H^n) at alpha = (1,2,5), eps = (3,7). For R the defect is 0 at (g,n) = (0,4), (1,1), (2,0) through y^2 and at (3,0) through y^1. For R^* the defect is non-zero already at y^1: 1249/2268, 76487/17280, -28742463168731/35831808000, -406623543620925998705062409/62912004762894336000000. Evidence only.""")

d.save()
