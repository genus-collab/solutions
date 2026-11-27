import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import Doc
d = Doc(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "LP2_REGULARITY_PROOF_draft.md"))
r = d.rep

OPEN_RULING = "Stage A uses eq:J=I and the genus-0 QDE at generic α on T_max, not on T_F. Ruling Q1 names statements on T_F. eq:J=I at T_max is a genus-0 mirror theorem, with no TR, cited by BS to Coates, Commun. Math. Phys. 287 (2009) 1071-1108 (MR2486673), not opened here. This use is outside the literal wording of Q1 and is recorded here pending the user's confirmation."

# header comment
r("""FRAMEB_PROOF.md as assembled and these edits are not refereed at the time of this edit: see s.0.4. -->""",
  """FRAMEB_PROOF.md as assembled and these edits are not refereed at the time of this edit: see s.0.4. Later in round 6c: minor items of the round-6c referee reports applied in place, where they could be re-derived. -->""")

# Rules
r("""Allowed: BS genus-0 mirror and oscillating-integral statements, cfasym, finite generation and algebraic independence; KP/Toda; topological recursion (TR) on the Hodge side. Forbidden: TR for the target theory, including the local-P2 mirror curve.""",
  """Allowed: BS genus-0 mirror and oscillating-integral statements; cfasym as a source of F2 only (not used, s.8.2); finite generation and algebraic independence; KP/Toda; topological recursion (TR) on the Hodge side. Forbidden: TR for the target theory, including the local-P2 mirror curve. Not inputs for the gap: BS rmk:refmm, BKMP/FLZ remodelling, HK10, Coates-Iritani 2018, BFGW Conj 1.9, BS prop:cfgraphsum, thm:cfasym, lem:cfpole.""")

# Prop 13.2_P2 label
r("""[P mod Stages A, B of FRAMEB_PROOF.md and Propositions 2.10, 2.12; inputs in s.8.2a]""",
  """[P mod Stages A, B of FRAMEB_PROOF.md, Propositions 2.10, 2.12 and PPZ Defs 2.2, 2.5, 2.13, Prop 2.14; for the clause F_cf = XX, the logarithmic form in Lemma 1.2(c) (DLMF 15.6.7, Gauss), not Barnes's second lemma; inputs in s.8.2a]""")

# Corollary G
r("""For g >= 2, with GW^cf_g := Cf(F_g^{T_F}),

    GW^cf_g = (-3)^{g-1} β_{2g-2}(ε1, ε2)/t_cf^{2g-2} + O(1).

O(1) means:""",
  """For g >= 2, with GW^cf_g := Cf(F_g^{T_F}), a germ at Δ_c = 0, and with t_cf the analytic continuation to a neighbourhood of Δ_c = 0 of BS's expression on (-1/27, 0) with log y = log(-y) + iπ, which equals t_c by Lemma 1.2(c),

    GW^cf_g = (-3)^{g-1} β_{2g-2}(ε1, ε2)/t_cf^{2g-2} + O(1).

O(1) means:""")

# G_{k,h}
r("""This is (Conifold asymptotics) of BS def:dirint. In particular F^cf_{0,h}""",
  """This is (Conifold asymptotics) of BS def:dirint. Here Π_{γ_cf} is read through its displayed expansion: BS's Π_{γ_cf} solves L2b∘θ = 0 and vanishes linearly at Δ_c = 0 with leading coefficient 1 (BS 929, 937, 955-957). Since L2b has exponents (0, 0) there, the solutions vanishing at Δ_c = 0 with no logarithm are the multiples of t_c, so Π_{γ_cf} = t_c = t_cf. The coefficient 𝔟_{k,h} depends on this normalisation; the integrality wording of BS 925-927 is not used. F -> F_cf in def:dirint is ev_cf on the polynomial of Theorem FB_B(e), whose coefficients lie in Q[X^{±1}, (1-X)^{-1}] rather than BS's Q[X^{±1}] (s.0.4). In particular F^cf_{0,h}""")

# referee status
r("""Each link had at least one referee pass in report form (s.0.4).""",
  """Each link had at least one referee pass, in report form or, for SP1-SP5 and Lemmas 2.1, 2.3-2.5, in an earlier round (s.0.4; FRAMEB_PROOF.md s.1.3).""")

# s.0.3
r("""and the interface statements Proposition 2.6, Theorem 2.7, Proposition 2.11. The labels (C1)-(C5)""",
  """and the interface statements Proposition 2.6, Theorem 2.7, Proposition 2.11. """ + OPEN_RULING + """ The labels (C1)-(C5)""")
r("""eq:J=I is cited by BS to MR2486673, which was not opened; it is not re-proved.""",
  """eq:J=I is cited by BS to MR2486673 (Coates, Commun. Math. Phys. 287 (2009) 1071-1108), which was not opened; it is not re-proved.""")
r("""  - Teleman, arXiv:0712.0160: Theorem (2), Prop 7.7, Prop 8.7(ii) (FB:A1.8, A2.8). The proof of Prop 7.7 was read, not re-derived; it is the one non-classical citation on which the graph sum rests.""",
  """  - Teleman, arXiv:0712.0160: Theorem (2), Prop 8.7(ii) (FB:A1.8), and 2.18, 2.19, that is Harer stability with Looijenga, and Madsen-Weiss (FB:A2.8). Since round 6c Prop 7.7 is not an input: FB:A2.8 derives the flatness equation directly, and (7.7.a) appears there only as a comparison. Theorem (2) is a non-classical citation on which the graph sum rests; eq:J=I is a further citation.""")
r("""  - Graber-Pandharipande (GP), arXiv:alg-geom/9708001: formula (1), Apps. A, C (FB:A1.2, B5).""",
  """  - Graber-Pandharipande (GP), arXiv:alg-geom/9708001: formula (1), Apps. A, C (FB:A1.2, B5). The fundamental-class and splitting axioms of the virtual class, with Coates-Givental App. 1 (iv)-(vi) for the twist (FB:A1.4, A1.5; Behrend-Manin not opened).""")
r("""The assembled file had its own referee pass in round 6 (L_logic: sound with minor fixes, no major item); see s.7.4 and s.8.2.""",
  """The assembled file had its own referee passes in round 6 (L_logic, and L_comp, independent computation: both sound with minor fixes, no major item); see s.7.4 and s.8.2.""")

# s.0.4
r("""(MV), Teleman's Prop 7.7 and eq:J=I stay cited, and their proofs were not re-derived.""",
  """(MV), Teleman's Theorem (2) with 2.18, 2.19, and eq:J=I stay cited, and their proofs were not re-derived.""")
r("""s.7.8 (D_VI); LEMMA_A_PROOF.md (L_logic).""", """s.7.8 (D_VI); LEMMA_A_PROOF.md (L_logic, L_comp).""")

# Lemma 1.2 remark
r("""in report form it had the passes refFB_C, refFB_C2. (c) is used only for Corollary G. The normalisation""",
  """in report form it had the passes refFB_C, refFB_C2. The logarithmic form in (c) is used for F_cf = XX in Proposition 13.2_P2; t_c = t_cf is used only for Corollary G. The normalisation""")

# Cor 2.2
r("""Thus PF fixes the row only up to a factor depending on z alone. (Round-3 fix: the factor is z b_1(z), not b_1(z).)""",
  """Thus PF fixes the row only up to a factor depending on z alone. (Round-3 fix: the factor is z b_1(z), not b_1(z).) If the b_n are constants and 𝒫(Sum b_n e_n) = 0 holds at every z-order (each order is a rational function of V), the same conclusion follows. By Theorem A, the expression is -(z^3/27) Sum_m C_m e_m with C_m := 27(m+1) b_{m+2} - 3(3m^2+3m+1) b_{m+1} + m^3 b_m. The terms with m < 0 vanish, since b_n = 0 for n <= 0. And e_m = (Vz/E)^m/prod_{j<=m}(1 + jz/E), so the coefficient of V^m is -(z^3/27) C_m (z/E)^m/prod_{j<=m}(1 + jz/E). Hence C_m = 0 for all m >= 0.""")

# s.2.2
r("""Not used in the proofs of Theorems LC_P2, R_P2. The items were called F1-F5""",
  """Not used in the proofs of Theorems LC_P2, R_P2; load-bearing for Proposition 13.2_P2 through FB:A3.5, B2, B4. The items were called F1-F5""")
r("""normalised formal Laplace expansion (BS line ~4245) at""",
  """normalised formal Laplace expansion (BS line 4256; Laplace expansion 4245-4252) at""")
r("""Result: d_p[e^{-φ_0/z} Λ] = e^{-φ_0/z} Λ(∂_p γ - z^{-1}(∂_p φ)γ). SP3 holds for any analytic amplitude γ.""",
  """Result: d_p[e^{-φ_0/z} Λ] = e^{-φ_0/z} Λ(∂_p γ - z^{-1}(∂_p φ)γ). Here the identity means (d_p - z^{-1} ∂_p φ_0) Λ_p(γ_p) = Λ_p(∂_p γ - z^{-1} (∂_p φ) γ). Likewise, in Lemma 2.4, L_pf(e^{u/z} f) stands for e^{u/z} 𝒫 f, which needs only θu = Y. SP3 holds for any analytic amplitude γ.""")
r("""Writing log F = Sum g_k w^k, (S) makes the difference of two solutions constant in β, and (D) at index k+1 kills those constants, so the solution is unique.""",
  """Each z-coefficient of F is a polynomial in β, because (1 + ξ/α)^{β-1} has jets polynomial in β; this is what makes ∂_β termwise legitimate in (D), where the term log α cancels. Write log F = Sum_{k>=1} g_k(β) w^k. The difference d_k of two polynomial solutions has period 1 by (S), so it is constant. (D) reads ∂_β g_1 = β - 1/2 and ∂_β g_{k+1} = -k g_k, which forces d_k = 0 for all k >= 1. For β not an integer, α^β needs a chosen log α; the applications have β in Z.""")

# Lemma 2.3
r("""For 𝒲, on a simply connected neighbourhood U of the semisimple locus in (y, α), the z^n coefficient of 𝒜 := Asym_{x^cr}[e^{𝒲/z} dΩ] lies in W^{-1/2} Q[α, Y][W^{-1}]. *Proof.* The z^n coefficient is the finite sum over m - k = n, 3k <= 2m, of (2^m m! k!)^{-1} (H^{ab}∂_a∂_b)^m [(-f_ng)^k g]|_{x_c}, polynomial in jets and in adj(H)/det H. The jets of 𝒲 of order |β| = m >= 2 are q_j δ_{β = m e_j} - Y/3^{m-1}. ∎""",
  """For 𝒲: let U be a simply connected open set of (y, α), y != 0, on which a root Y of L^ is analytic with Y W(Y) prod_j (Y + α_j) != 0. On U, 𝒜 := Asym_{x^cr}[e^{𝒲/z} dΩ] = c_0 W^{-1/2} Φ(α, Y; z), with c_0^2 = -1/3 and [z^n] Φ in Q[α, Y, W^{-1}], Φ_0 = 1. *Proof.* Write φ = -𝒲, with H its Hessian, and γ for the amplitude. The z^n coefficient is the finite sum over m - k = n, 3k <= 2m, of (2^m m! k!)^{-1} (H^{ab}∂_a∂_b)^m [(-φ_ng)^k γ]|_{x_c}, polynomial in jets and in adj(H)/det H. The jets of 𝒲 of order |β| = r >= 2 are q_j δ_{β = r e_j} - Y/3^{r-1}. The constant c_0 comes from dΩ = (1/3) d^3x and det Hess 𝒲 = W/3 (Lemma 2.1). ∎""")

# Lemma 2.5
r("""*Proof.* (i): (Y + zθ)^n F = Y^n F + z[n Y^{n-1} θF + C(n,2) Y^{n-2} (θY) F] + O(z^2); use L^' = -W/Y and θY = -Y prod(Y+α)/W. (ii):""",
  """*Proof.* (i): put T := zθ + Y. Then 𝒫 = prod_j (T + α_j) + 27 y T^3 + y(27 z T^2 + 6 z^2 T), because 3T(3T+z)(3T+2z) = 27T^3 + 27zT^2 + 6z^2T. Since T^n F = Y^n F + z[n Y^{n-1} θF + C(n,2) Y^{n-2} (θY) F] + O(z^2), we get 𝒫_0 = L^(Y) = 0 and 𝒫_1 F = L^'(Y) θF + [½ L^''(Y) θY + 27 y Y^2] F, where ' is the Y-derivative at fixed y. Now L^' = -W/Y and θY = -Y prod(Y+α_j)/W = 27 y Y^4/W. Differentiating W = -Y L^'(Y) along θ gives θW = (W/Y) θY - Y L^'' θY - 81 y Y^3. These give ½ L^'' θY + 27 y Y^2 = -θW/(2Y), which is (i). (ii):""")
r("""(iii): evaluate. ∎

**Proposition 2.6""",
  """(iii): evaluate. ∎

The proof is algebraic. (i) and (ii) hold verbatim with O(Ω') replaced by F[[y]], with F a field of characteristic 0, and Y in F[[y]] with L^(Y) = 0 and Y(0) W(Y)(0) != 0. The θ-constants are then F, and f, g are series in z with coefficients in F[[y]].

**Proposition 2.6""")

# Prop 2.6
r("""uses SP1-SP5. (C2) is discharged: FB:A2.5, A3.3.] Fix generic α""",
  """uses SP1-SP5. (C2): as a statement about the row eq:canR, (Normalisation) is its definition, and this label is complete. That this row is Teleman's R adds the inputs of FB:A2.8, A2.9, A3.2, A3.3 (Teleman, PPZ, Mumford as in FP, FP Theorem 4).] Fix generic α""")
r("""where c in C^× is a universal constant with c^2 = ±1, fixed by the Gaussian square-root convention.""",
  """where c in C^× is a universal constant, with c^2 = -1 for BS's Asym, normalised by (det Hess(-𝒲))^{-1/2}, and c = 1 if (det Hess 𝒲)^{-1/2} is used instead.""")

# Thm 2.7
r("""the reduced (C1) of s.0.3; (C2), (C3) are discharged: FB:A2.5, A2.6, A3.3; uses SP1-SP5. For the role""",
  """the reduced (C1) of s.0.3; (C3) is discharged: FB:A2.6; (C2): as a statement about the row eq:canR, (Normalisation) is its definition, and this label is complete. That this row is Teleman's R adds the inputs of FB:A2.8, A2.9, A3.2, A3.3 (Teleman, PPZ, Mumford as in FP, FP Theorem 4). Uses SP1-SP5. For the role""")

# Scope s.2.7
r("""still needs Theorem 2.7(e), now in the form FB:B2, B4, hence SP1-SP5""",
  """needs FB:B2, B4 (the content of Theorem 2.7(e)), hence SP1-SP5""")

# Lemma A refereeing, line 822
r("""The assembled LEMMA_A_PROOF.md had its own referee pass in round 6 (L_logic: sound with minor fixes, no major item), which answers its open point O1.""",
  """The assembled LEMMA_A_PROOF.md had its own referee passes in round 6 (L_logic, and L_comp, independent computation: both sound with minor fixes, no major item), which answer its open point O1.""")

# B-c
r("""but the Teleman input of conifold Prop 13.2 was read through summaries only (conifold critic gap list), and it keeps a BS citation.""",
  """but the Teleman input of conifold Prop 13.2 stays a citation (arXiv:0712.0160v3 read in round 6), its field-of-definition paragraph carries an OPEN mark from round 6c (LSP Prop 13.2(a), sign of u), and it keeps a BS citation.""")

# Cor G proof
r("""For the Taylor coefficients: wc, t_c and F_cf lie in Q[[Δ_c]] and y^{-1} is regular at Δ_c = 0, so an element of Q(ε1, ε2)[y^{±1}, Δ_c^{-1}][F_cf] has a Laurent expansion in Δ_c over Q(ε1, ε2), and t_c is a coordinate over Q. The difference is such an element, up to the constant, and it has no negative powers. ∎""",
  """For the Taylor coefficients: wc, t_c and F_cf lie in Q[[Δ_c]] and y^{-1} is regular at Δ_c = 0, so an element of Q(ε1, ε2)[y^{±1}, Δ_c^{-1}][F_cf] has a Laurent expansion in Δ_c over Q(ε1, ε2), and t_c is a coordinate over Q. GW^cf_g is such an element, up to the constant, and t_c^{2-2g} lies in Q((Δ_c)). So the difference has a Laurent expansion in Δ_c over Q(ε1, ε2), equivalently in t_c, and by Theorem R_P2 and Corollary SV_P2 it has no negative powers. F_cf and t_c converge near Δ_c = 0, t_c is a local coordinate, and y^{-1} is regular there. So the difference is a finite expression that is meromorphic on a disc about Δ_c = 0, and with no negative powers it is holomorphic. ∎""")
r("""and uniqueness of the polynomial in F gives Cf(GW_g) = Sum_j E^{2g-2j} (ε1 ε2)^{j-1} Cf(f_j).""",
  """and uniqueness of the polynomial in F gives Cf(F_g^{T_F}) = Sum_j E^{2g-2j} (ε1 ε2)^{j-1} Cf(f_j); with BS's d > 0 convention the left side changes by the constant GW^{T_F}_{g,0}, which affects only the O(1) term.""")

# s.8.1 table
r("""    Lemma 1.2(c) (t_c = t_cf, A_cf = C1cf, branch)   P mod C     DLMF 15.6.7, Barnes's second lemma, Gauss's ψ values; FB:C Step 6;
                                                                 used only for Corollary G""",
  """    Lemma 1.2(c) (t_c = t_cf, A_cf = C1cf, branch)   P mod C     DLMF 15.6.7, Barnes's second lemma, Gauss's ψ values; FB:C Step 6;
                                                                 its logarithmic form is used for F_cf = XX in Proposition 13.2_P2;
                                                                 t_c = t_cf is used only for Corollary G""")
ROUTE = "not used in Theorems LC_P2, R_P2;\n" + " "*65 + "used for Proposition 13.2_P2 and Corollary G (FB:A2.5, A3.5, B1, B2, B4)"
r("""S (written) R3 Part C; interface only""", """S (written) R3 Part C; """ + ROUTE)
r("""P           R2:Lemma 1; interface only""", """P           R2:Lemma 1; """ + ROUTE)
r("""P           R2:Lemma 2; interface only""", """P           R2:Lemma 2; """ + ROUTE)
r("""P           SP2, SP3; interface only""", """P           SP2, SP3; """ + ROUTE)
r("""P           R2:Lemma 4; interface only""", """P           R2:Lemma 4; """ + ROUTE)
r("""P mod C     eq:J=I, genus-0 QDE; SP1-SP5; (C2) discharged (FB:A2.5, A3.3);
                                                                 R2:Prop 5; formal form FB:A3.5; interface only""",
  """P mod C     eq:J=I, genus-0 QDE; SP1-SP5; (C2): label in s.2.4 (that the row is Teleman's R adds
                                                                 the inputs of FB:A2.8, A2.9, A3.2, A3.3); R2:Prop 5; formal form FB:A3.5;
                                                                 not on the route of Proposition 13.2_P2 (FB:B2, B4, B6 replace it); kept as the analytic statement""")
r("""P mod C     eq:J=I, genus-0 QDE; SP1-SP5; (C2), (C3) discharged (FB:A2.6);
                                                                 R2:Thm 6, Prop 8; R3:Cor B; steps 2, 5 of (e) as in FB:B4; interface only;
""",
  """P mod C     eq:J=I, genus-0 QDE; SP1-SP5; (C3) discharged (FB:A2.6); (C2): label in s.2.4;
                                                                 R2:Thm 6, Prop 8; R3:Cor B; steps 2, 5 of (e) as in FB:B4;
                                                                 not on the route of Proposition 13.2_P2 (FB:B2, B4, B6 replace it); kept as the analytic statement;
""")
r("""answered by the round-6 referee L_logic (no major item);""",
  """answered by the round-6 referee L_logic, with L_comp (independent computation) as second pass (no major item);""")
r("""P mod (MV)  LLZ eq. (1) only, at μ = 1^a; LEMMA_A_PROOF.md Prop 2.1, Lemma 3.2, (4.2);""",
  """P mod (MV)  LLZ eq. (1) only, at μ = 1^a; LEMMA_A_PROOF.md Prop 2.1, Lemma 3.2, (4.2) at k = 2;""")
r("""P mod C     eq:J=I, genus-0 QDE; Teleman Thm (2), Prop 7.7, Prop 8.7(ii); PPZ; GP; Mumford, FP Thm 4;
                                                                 A3.5 also SP1-SP5, Lemmas 2.1, 2.3-2.5; two referee passes on the report form""",
  """P mod C     eq:J=I, genus-0 QDE, eq:SHb; Teleman Thm (2), Prop 8.7(ii), 2.18, 2.19; PPZ; GP; virtual-class axioms; Mumford, FP Thm 4;
                                                                 A3.5 also SP1-SP5, Lemmas 2.1, 2.3-2.5; two referee passes on the report form;
                                                                 A2.8 as rewritten in round 6c (no Prop 7.7): no further pass""")
r("""P mod C     Stages A, B; Props 2.10, 2.12; FB:C Step 5; s.8.2a""",
  """P mod C     Stages A, B; Props 2.10, 2.12; PPZ; logarithmic form of Lemma 1.2(c) for F_cf = XX; FB:C Step 5; s.8.2a""")

# s.8.2a
r("""Dependency: Stage A -> Stage B -> Stage C (Proposition 13.2_P2). Lemma (t_cf) is independent of Stages A, B. Corollary G = Proposition 13.2_P2 + Lemma 1.2(c) + Theorem R_P2 + Corollary SV_P2 + Theorem VI.""",
  """Dependency: Stage A -> Stage B -> Stage C (Proposition 13.2_P2). Lemma (t_cf) is independent of Stages A, B. Its logarithmic form gives F_cf = XX in Proposition 13.2_P2; its constant (Barnes) gives t_c = t_cf in Corollary G. Corollary G = Proposition 13.2_P2 + Lemma 1.2(c) + Corollary 2.8 + Theorem R_P2 + Corollary SV_P2 + Theorem VI.""")
r("""eq:J=I is cited by BS to MR2486673, not opened. Definitions only:""",
  """eq:J=I is cited by BS to MR2486673 (Coates, Commun. Math. Phys. 287 (2009) 1071-1108), not opened. Definitions only:""")
r("""Asym at about 4245, t_cf at 5145.""", """Asym at 4256 (Laplace expansion 4245-4252), t_cf at 5145.""")
r("""2. **Teleman**, Theorem (2), Prop 7.7, Prop 8.7(ii). **PPZ**, Defs 2.2, 2.5, 2.13, Prop 2.14. **GP**, formula (1), Apps. A, C. **Mumford**""",
  """2. **Teleman**, Theorem (2), Prop 8.7(ii), 2.18, 2.19 (Harer stability with Looijenga; Madsen-Weiss). **PPZ**, Defs 2.2, 2.5, 2.13, Prop 2.14. **GP**, formula (1), Apps. A, C. The fundamental-class and splitting axioms of the virtual class, with Coates-Givental App. 1 (iv)-(vi) for the twist (Behrend-Manin not opened). **Mumford**""")
r("""Gauss's values of ψ at 1/3, 2/3.** Used only for Lemma 1.2(c).""",
  """Gauss's values of ψ at 1/3, 2/3.** Used for Lemma 1.2(c), which enters Proposition 13.2_P2 (logarithmic form only) and Corollaries G, G_{k,h}.""")
r("""What the referee passes did not cover is FRAMEB_PROOF.md s.1.4: the proofs of Teleman's Prop 7.7 and Theorem (2) (read only);""",
  """What the round-6b passes did not cover is FRAMEB_PROOF.md s.1.4: the proofs of Teleman's Theorem (2) and 2.17-2.19 (used as stated);""")

# s.8.3
r("""with status [P mod Stages A, B of FRAMEB_PROOF.md and Propositions 2.10, 2.12]. Its proof""",
  """with the status given there. Its proof""")

# errata
r("""12. **Lines 4055-4057.** Regularity""", """12. **Lines 4054-4056.** Regularity""")
r("""and item 9 (lines 4071-4076: regularity at y = 0 is argued for the coefficients of the polynomial in F, while what is given is regularity of the total, and F = O(y); no counterexample is known, s.0.4).""",
  """, item 9 (lines 4071-4076: regularity at y = 0 is argued for the coefficients of the polynomial in F, while what is given is regularity of the total, and F = O(y); no counterexample is known, s.0.4) and item 11 (lines 925-927: the expansions of Π_γ are normalisation statements, and hold for an integral cycle only up to a constant prefactor).""")

# Appendix
r("""Copies of the round 1-3 scripts are in LP/round1_scripts, LP/round2_scripts, LP/round3_scripts.""",
  """Copies of the round 1-5 scripts are in LP/round1_scripts ... LP/round5_scripts, and of round 6 in LP/round6_scripts/<key>.""")
r("""**Round 6 (SP/lp2_round6/: D_rows, D_core, D_hodge, D_VI, L_logic).**""",
  """**Round 6 (SP/lp2_round6/: D_rows, D_core, D_hodge, D_VI, L_logic, L_comp).**""")

d.save()
