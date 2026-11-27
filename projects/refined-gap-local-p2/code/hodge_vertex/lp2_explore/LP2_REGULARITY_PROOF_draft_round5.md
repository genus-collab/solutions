<!-- Round 5 (2026-11-02): edited in place from the round-4 consolidated draft (backup: LP2_REGULARITY_PROOF_draft_round4.md). Applies ROUND4_CRITICS_2026-11-02.md item by item, integrates the round-5 reports "modrel" and "VI" with their referee fixes, and points Lemma A to hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md. Built on ROUND2_PROOFS_rows_on_TF.md, ROUND3_PROOFS_H2_edges_PF.md and the status notes. Journals (under <user home>/.claude/projects/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/subagents/workflows/<run>/journal.jsonl): round 1 wf_e928f771-9e3 (this is the only record of the round-1 frame-B skeleton and of the items N1-N8 quoted below, except N4), round 2 wf_1ab43e48-d73, round 3 wf_74ed0422-50f, round 4 wf_bab37015-382 (reports thmR, lemmaA and their referees), round 5 wf_9d7fd94e-8d7. NOT yet refereed as a whole. -->
# Regularity at the conifold of the frame-B Hodge-vertex graph sum for local P2 x C^2, to all orders and in all genera

K_P2 x C^2 on the fibrewise torus T_F, conifold point y = -1/27, conifold polarisation, Hodge-vertex frame B, mode B.
Modelled on hodge_vertex_conifold/leaf_proof/leaf_structure_proof.md (cited as LSP).

**Labels.** Every statement carries one of the following labels.
- **[P]** proved in this text, to all orders. Where a proof is taken over from an earlier round, the round and item are given (R2:x = ROUND2_PROOFS_rows_on_TF.md item x; R3:x = ROUND3_PROOFS_H2_edges_PF.md item x; the round-4 reports thmR, lemmaA and their referees are in ROUND4_CRITICS_2026-11-02.md and journal wf_bab37015-382; the round-5 reports modrel, VI, ThmL, MVmatch are in journal wf_9d7fd94e-8d7).
- **[C]** cited. The citation gives the file and section, or the BS v2 label with its line in BS_v2_2026-11-02/main.tex.
- **[S]** standard and written out here with its proof: formal stationary phase (SP1-SP5), s.2.2. An [S] item is not an open input.
- **[O]** open.
- **[E]** finite-order computer evidence. It is never used as a proof step (Appendix).

**Rules (user rulings, ALLGENUS_PATHWAY_2026-11-01.md s.6).** Allowed: BS genus-0 mirror and oscillating-integral statements, cfasym, finite generation and algebraic independence; KP/Toda; topological recursion (TR) on the Hodge side. Forbidden: TR for the target theory, including the local-P2 mirror curve. Everything must hold to all orders: finite-order computations are evidence only and may not fix a formula or normalisation that a proof relies on. The gap corollary is in scope as a final step. No step below uses TR for the target theory. Hodge-side TR enters only through the backup route for Lemma A (s.7.4) and the cross-check of Theorem VI (s.7.8), and neither is a proof step.

Contents. 0 Summary. 1 Conventions. 2 Rows on T_F, including s.2.7 (lem:modrel and the interface with BS). 3 The singular colour: F1, F2, L1, Z-infinity, closure, K = 1, Dil_0. 4 Proposition H2. 5 The residue lemma and E_00. 6 Regular colours and mixed edges (Theorem M). 7 Hodge side and Theorem R_P2, including s.7.8 (proof of the vertex identity). 8 Status of every statement; remaining inputs. 9 Errata in BS v2. Appendix: computer checks.


## 0. Summary

### 0.1 The objects

Fix ε1, ε2 with ε1 ε2 E != 0, where E := ε1 + ε2. In the code these are called e1, e2. The Hodge weights are w = (E, -ε1, -ε2). For g >= 2, the frame-B Hodge-vertex graph sum of genus g is

    G_g := Sum_Gamma |Aut Gamma|^-1 Sum_{c: V -> {0,1,2}} T(Gamma, c),
    T(Gamma, c) := Sum_dec prod_v [ Delta_{c_v}^{g_v - 1} (1/m_v!) Int_{Mbar_{g_v, n_v + m_v}} Omega_w prod_{h at v} ψ_h^{a_h} prod_{k=1}^{m_v} ψ_{n_v+k}^{d_{v,k}} · prod_k Dil^_{c_v}(d_{v,k}) ] prod_e Ehat_{c c'}(a_e, b_e).      (0.1)

The ingredients are as follows.
- Gamma runs over connected stable graphs of genus g without legs; n_v is the number of half-edges at v.
- **Decorations; (0.1) is their definition.** A decoration assigns an exponent a_h >= 0 to each half-edge (a_e, b_e denote the exponents at the two half-edges of e), and to each vertex an integer m_v >= 0 together with an ordered tuple (d_{v,1}, ..., d_{v,m_v}) of dilaton degrees >= 2, with weight 1/m_v!. Since Omega_w has mixed degree, the integral vanishes unless Sum a_h + Sum d_{v,k} <= dim = 3g_v - 3 + n_v + m_v. As d >= 2 this forces m_v <= 3g_v - 3 + n_v. So the decoration sum is finite, and T(Gamma, c) is a finite sum of products of the factors below; this is what gives "analytic at the conifold" its meaning in Theorem R_P2.
- Omega_w is the normalised Hodge class prod_{a=1}^{3} Lambda^vee_{g_v}(w_a)/w_a^{g_v} (s.1.5). It is the same at every colour and independent of y.
- Delta_i, Dil^_i and Ehat_ij are the normalised frame-B factors of s.1.4. They are built from the rows of Definition 1.3, taken in the conifold polarisation c1 = C1cf.

Colour 0 is the singular colour: on the real branch t = y^{1/3} -> -1/3 it is the unique colour whose critical point escapes (Lemma 2.1).

(0.1) is the definition of the decoration sum. It reproduces the code's value (Lemma 1.4). Theorems LC_P2 and R_P2 are statements about this sum. Identifying G_g with a Gromov-Witten potential is a separate matter (Proposition 13.2_P2, s.8.3).

### 0.2 The main theorems

**Theorem LC_P2 (leaf structure at the singular colour; F1-F3, L1-L3, to all orders in t_c and z).** [P: no BS input, since Proposition 2.10 is now proved from the closed form]

Hypotheses: ε1, ε2 complex with ε1 ε2 E != 0, or indeterminates. All factors are those of s.1.4 at c1 = C1cf. Notation: ρ := -E v_0/t_c (s.1.2); λ_j(ψ) := 1/(E - jψ); Λ_m := span{λ_1, ..., λ_m}; q_t(ψ) := Sum_{d>=1} Dil^_0(d) ψ^d; q_0 := E - E^2 λ_1. Then:

- (F1) Delta_0 = 3 ε1 ε2 E^3/(t_c ρ)^2. Here ρ is an analytic unit at t_c = 0 with ρ(0) = 1. For j = 1, 2, Delta_j is analytic at t_c = 0, with Delta_j(0) = ε1 ε2 E^3 ω^{2j}/(3(1-ω^j)^2) != 0.
- (F2) The normalised identity row at colour 0 satisfies P^_0^{(0)}(z; t_c = 0) = E/(E+z). The rescaled row σ := ρ P^^{(0)}(-ψ) has σ^(-1) = (0, E, -E^2) and σ^(0) = (E λ_1, (E^2/9) λ_1, 0).
- (F3) Every dilaton factor Dil^_i(d) (i in {0,1,2}, d >= 1) and every edge factor Ehat_ij(m, n) (i, j in {0,1,2}, m, n >= 0) is analytic at t_c = 0.
- (L1) [t_c^0](ρ q_t) = q_0. For every k >= 1, [t_c^k](ρ q_t) lies in Λ_{k+1}; in particular it has no constant term.
- (L2) Let j in {1, 2} and let e, b, k be integers with b, k >= 0. Then [t_c^k][w^b](ρ^e Ehat_0j(ψ, w)) lies in Λ_{k+1}(ψ) and has no constant. Moreover Ehat_j0(w, ψ) = Ehat_0j(ψ, w), so this covers both orientations of a mixed edge. For an edge whose other end is at colour 0: [t_c^k][ψ'^b](ρ^e Ehat_00) lies in Λ_{k+1}(ψ).
- (L3) For every integer e and every k >= 0, [t_c^k](ρ^e Ehat_00(ψ, ψ')) lies in span{λ_i(ψ) λ_l(ψ') : i, l >= 1, (i-1)+(l-1) <= k} and is symmetric. Furthermore Ehat_00|_{t_c=0} = -(ε1 ε2 E^4/27) λ_1 λ_1'.

Every membership is an identity of formal power series in ψ (and ψ', w), holding at every ψ-order and every t_c-order. The coefficients of the leaves at each order lie in Q(ε1, ε2), up to the algebraic factors in ω at the regular colours.

**Theorem R_P2 (regularity of graphs with an edge).** [P mod Lemma A(h), 2 <= h <= g-1]

Fix g >= 2. Assume M_irr(h) (Lemma A; s.1.5 and s.7.4) for 2 <= h <= g-1. Then for every connected stable graph Gamma of genus g with at least one edge, and for every colouring c, the decoration sum T(Gamma, c) of (0.1) is analytic at y = -1/27, i.e. its Laurent expansion in t_c has no negative powers. The same therefore holds for the colour sum.

**Corollary SV_P2 (the single-vertex graph).** [P mod Lemma A(h) for 2 <= h <= g, and Theorem VI, which is P mod LLZ eqs. (1), (5) and the boson-fermion correspondence (s.7.8)]

Fix g >= 2 and assume M_irr(h) for 2 <= h <= g. The single-vertex graph is analytic at colours 1 and 2. At colour 0 its polar part is the single term

    (3 ε1 ε2 E^3)^{g-1} Z[g,0] t_c^{2-2g},   Z[g,0] := Fpot_g(q_0).

By the vertex identity (VI), Theorem VI of s.7.6, this equals (-3)^{g-1} β_{2g-2}(ε1, ε2) t_c^{2-2g}, where β_n := (n-1)! [t^n] 1/(4 sinh(ε1 t/2) sinh(ε2 t/2)) (BS eq. at line 959). Hence the polar part of G_g is exactly (-3)^{g-1} β_{2g-2} t_c^{2-2g}.

**Corollary G (refined gap, BS rmk:cfgap, line 5422).** [P mod the above and Proposition 13.2_P2 [O]]

If in addition G_g = GW^cf_g (Proposition 13.2_P2, open), then GW^cf_g = (-3)^{g-1} β_{2g-2}/t_cf^{2g-2} + O(1). The identification t_c = t_cf is Lemma 1.2(c). This would prove BS's conjecture rmk:cfgap. **This corollary is NOT claimed** while Proposition 13.2_P2 is open.

### 0.3 What is cited
- **BS genus-0 inputs (allowed).** Since round 5 these are **not** inputs of Theorems LC_P2 and R_P2, or of Corollary SV_P2 on the recorded route for (VI): Proposition 2.10 is proved from Theorem A alone, and s.3-s.6 contain no other BS link (dependency reading by the round-5 modrel author and referee). (η^{ab} and L_pf enter those theorems as definitions.) The BS inputs are needed only to identify the rows and the sum with BS's objects, i.e. for Lemma 1.2(c), Proposition 2.6, Theorem 2.7, Proposition 2.11, Corollary G and Proposition 13.2_P2. They are labelled (C1)-(C3) as in R2:s.0; the parenthesised labels (Cn) are never to be read as C1cf or C1i.
  - (C1) eq:Ifun 2621, eq:J=I 2663, eq:asymp 2943 (with about 2950) and eq:canR 3566: for generic α, e^{u_i/z} R^α_{0i} formally solves eq:pf (2628).
  - (C2) (Normalisation) 3497, eq:phi 3031 and the weights at about 2383: the value at y = 0 and Psi_{0i} = Delta_i^{-1/2}.
  - (C3) eq:RHb 3573, eq:ADEFXdef 2681, eq:etaflat 2419 and (Symplecticity) 3496. lem:modrel 2695 is no longer cited: both of its relations are proved in Proposition 2.12.
  - Also used: (C4) prop:algindABC 2854, only to make Cf well defined; and t_cf at line 5145 with eq:pfinf 2757, only for Lemma 1.2(c).
  - For Proposition 2.11 and Proposition 13.2_P2 only: (C5) eq:gtH 3799, thm:haegen 3902, lem:lho 3247 and the regularity on the diagonal at about line 4056.
- **Lemma A (Hodge side).** This is the conifold Lemma A, transferred in s.7.4. Its written proof is hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md: Theorem L for every h >= 2 and every non-empty R, proved there modulo the Marino-Vafa formula (LLZ eq. (1); OP Theorem 1 is a second published proof) and textbook facts. The assembled file has not yet had its own referee; see s.7.4 and s.8.2.
- **Hodge-side standard facts.** These are the string and dilaton equations for Omega_w; they hold because Lambda^vee is pulled back under forgetful maps.
- **Vertex identity (VI).** Used only for Corollary SV_P2. It is now proved in all genera (Theorem VI, s.7.8) from LLZ eqs. (1) and (5), the boson-fermion correspondence, the string and dilaton equations and (7.2). It no longer rests on the ABDKS/Eynard normalisation.

### 0.4 What is not claimed
- The identification of G_g with a GW potential. This is Proposition 13.2_P2, s.8.3.
- That the code (genRmat, cPP_z9) computes the rows of Definition 1.3 to all orders. That is evidence only (s.2.5, Appendix).
- Corollary G.
- A literature-free proof of Lemma A or of Theorem VI: both are proved modulo the literature listed in s.8.2.
- That the whole text has been refereed. Proposition 2.10, Proposition 2.12 and Theorem VI have had one referee pass each (sound with minor fixes, applied here); the extension of Fpot_h in s.1.5 and (H0) of s.7.1 are new in this edit and unrefereed.


## 1. Conventions

### 1.1 Parameters and variables
- **Weights.** ε1, ε2 (the code's e1, e2), E := ε1 + ε2, with ε1 ε2 E != 0. The elementary symmetric functions of α = (α0, α1, α2) are written e_k(α), never e_k alone. T_F is the point α* = (E/3)(1,1,1). On T_F we write α := E/3.
- **Base coordinate.** θ := y d/dy, t := y^{1/3}, ω := e^{2πi/3}, Δ_c := 1 + 27y. Note θ = (Δ_c - 1) d/dΔ_c.
- **Colours.** t_i := ω^i t for i = 0, 1, 2, and
  - Y_i := -α/(1+3t_i),
  - q_i := Y_i + α = E t_i/(1+3t_i),
  - v_i := 1/q_i = (1+3t_i)/(E t_i) = 3/(E+3Y_i),
  - V_i := E v_i.

  Caution: the script edges_sing.pl (line 129, variable $vi) uses v for 1 + 3t_i. That is a different quantity, and it is never used here.
- **Metric.** W(Y) := e_1(α)Y^2 + 2e_2(α)Y + 3e_3(α) (eq:delta1, 2971). On T_F: W(Y_i) = E q_i^2 = E^3 t_i^2/(1+3t_i)^2, Delta_i := 3 ε1 ε2 W(Y_i) = 3 ε1 ε2 E/v_i^2, and we fix Delta_i^{1/2} := sqrt(3 ε1 ε2 E)/v_i. This is the code's sqDsym, and it is transported together with the rows under t -> ω^i t.
- **Row and edge variables.** z is the argument of a row. ψ, ψ', w are edge variables: a half-edge with variable ψ carries the row at z = -ψ, i.e. P^^{(i)}(-ψ). So Ehat_ij(ψ, w) of (1.1) has edge variables as arguments, not row arguments.
- **Leaves.** λ_j(ψ) := 1/(E - jψ) for j in Z \ {0}, and λ_0 := 1/E. Further Λ_m := span{λ_1, ..., λ_m}, Λ_0 := 0, Λ_m^+ := C + Λ_m; in bounds, Λ_m is read as Λ_{max(m,0)}.
  **Normalisation warning.** R2:Cor 9 uses (1 - jψ/E)^{-1} = E λ_j, while R3 and this text use λ_j. Spans are unaffected, but the coefficient statements σ_0^(0) = E λ_1 and Sum_j β_j/j = E δ_{k0} are correct only in the present normalisation.
- **Leaf identities (all elementary).**
  - (LAM1) -ψ λ_j = (1 - E λ_j)/j, for j != 0.
  - (LAM2) (E - Nψ) λ_j = N/j + (1 - N/j) E λ_j, for j != 0; in particular (E - Nψ) λ_N = 1 for N != 0.
  - (LAM3) λ_l(ψ) - λ_l(-w) = l (ψ + w) λ_l(ψ) λ_l(-w).
  - (LAM4) E λ_1 = 1 + ψ λ_1.
  - The base point is q_0 := E - E^2 λ_1 = -Eψ/(E-ψ), and T := q_0 + ψ = -ψ^2/(E-ψ).

### 1.2 The conifold coordinate and the unit ρ
L2b := θ^2 + 3y(3θ+1)(3θ+2) = Δ_c θ^2 + 27yθ + 6y. Its exponents at Δ_c = 0 are (0, 0), because near Δ_c = 0 it reads Δ_c D^2 + D + ... with D = d/dΔ_c. Let wc be the solution analytic at Δ_c = 0 with wc(0) = 1, and let t_c be defined by θ t_c = -wc, t_c(Δ_c = 0) = 0. Put C1cf := θ t_c = -wc. For any section c1 put C1i := 1/c1, XX := θ log c1, d2 := E(C1i^2 - 1)/(2Δ_c). This agrees with lp2_lib.m tcSetup lines 191-201.

**Lemma 1.2.** [P]
- (a) t_c = Δ_c + O(Δ_c^2) is an analytic coordinate at Δ_c = 0, and C1cf(0) = -1. At c1 = C1cf, the functions C1i, XX and d2 are analytic at Δ_c = 0: d2 is analytic because C1i(0)^2 = 1.
- (b) On the real branch t = -(1 - Δ_c)^{1/3}/3 (colour 0): v_0 = -t_c/E + O(t_c^2). So ρ := -E v_0/t_c is an analytic unit with ρ(0) = 1, and ρ^2 Delta_0 = 3 ε1 ε2 E^3/t_c^2.
- (c) t_c = t_cf and C1cf = A_cf, with t_cf, A_cf as in BS sec:cfgen, line 5145.

*Proof.*
- (a) From the definition: dt_c/dΔ_c = -wc/(Δ_c - 1), which equals 1 at 0.
- (b) Δ_c = (1+3t)(1-3t+9t^2) and the second factor equals 3 at t = -1/3, so 1+3t = Δ_c/3 + O(Δ_c^2) and v_0 = (1+3t)/(Et) = -Δ_c/E + O(Δ_c^2). The identity follows from Delta_0 = 3 ε1 ε2 E/v_0^2.
- (c) t_cf is a combination of the periods of eq:pfinf (2757), whose operator is L2b∘θ, and BS 5145 states that t_cf = Δ_c + O(Δ_c^2), i.e. it is analytic and vanishes at Δ_c = 0. Then θ t_cf is an analytic solution of L2b. Since the exponents are (0,0), the analytic solutions form the line C·wc. So θ t_cf = κ θ t_c, and t_cf = κ t_c + const. Vanishing kills the constant, and the linear coefficient gives κ = 1. ∎

Remark: (c) was new in round 4. Critics 1 and 2 of round 4 read eq:pfinf (it is L2b∘θ) and BS 5145 and agree with the proof, but it has had no separate referee pass. It is used only for Corollary G. The normalisation of t_c enters every exact constant through ρ(0) = 1, and through Z-infinity (Corollary 2.9, s.3.1). All filtration statements hold for any ρ(0) != 0.

### 1.3 Rows
Special functions:
- Θ(x) := exp Sum_{m>=1} B_{2m} x^{2m-1}/(2m(2m-1)), so Θ(-x) = 1/Θ(x). The two-variable version used in SP5 and Theorem 2.7 is BS's (line 4278): Θ(z, x) := exp[-Sum_{k>=1} B_{k+1}(x) (-z)^k/(k(k+1))], with B_n(x) the Bernoulli polynomials. Then Θ(z, 0) = Θ(z), Θ(z, x+1) = (1 + xz) Θ(z, x), Θ(z, -m) = Θ(z)/prod_{j<=m} (1 - jz);
- h_w(z) := Θ(-z/E) Θ(z/ε1) Θ(z/ε2), with h_w(z) h_w(-z) = 1 (frame-B twist; this is Υ(z/ε1, z/ε2) of BS line 4282);
- N(z) := sqrt(ε1 ε2) Θ(-z/ε1) Θ(-z/ε2);
- the recurrence (Rec): r_{-1} = 0, r_0 = 1, 27(m+1) r_{m+1} = 3(3m^2+3m+1) r_m - m^3 r_{m-1} for m >= 0 (r = 1, 1/9, 2/81, 14/2187, ...; BS v2 displays only r_0, r_1, r_2; 14/2187 is computed);
- S(v, z) := Sum_{l>=0} r_l (vz)^l / prod_{j=1}^{l+1} (1 + jz/E);
- S^[m](v, z) := Sum_l r_l (vz)^l / prod_{j=1}^{l+1-m} (1 + jz/E), where a product with upper limit < 1 is read as 1 (the case l = 0, m = 2 occurs in Lemma 4.1).

Each z-coefficient of S is a polynomial in v, and [z^k]S has degree <= k in v. This is immediate from the definition, since only l <= k contribute.

**Definition 1.3 (the rows of the theorems).** At colour i and for a non-vanishing section c1(y):
- P_0^{B,(i)} := Delta_i^{-1/2} S(v_i, z);
- P_1^{B,(i)} := C1i (Y_i + zθ) P_0^{B,(i)};
- P_2^{B,(i)} := (Δ_c/C1i^2) [(Y_i + zθ) P_1^{B,(i)} - d2 P_1^{B,(i)}].

The large-radius-frame rows are P_a^{(i)} := h_w(z) P_a^{B,(i)}. The normalised rows are P^_a^{(i)} := Delta_i^{1/2} P_a^{B,(i)}, so that **P^_0^{(i)} = S(v_i, z) is component 0 of the colour-i normalised row, not the whole row.** The theorems take c1 = C1cf, the same section at every colour; this is legitimate because c1 depends on y only.

Two remarks.
- Definition 1.3 uses only XX = θ log c1 (θ C1i = -C1i XX), so L2b plays no role in defining the rows. L2b is needed only for Cf to commute with θ on the whole ring.
- With c1 = A these are eq:RHb (3573), through Δ_c A^2 = 1/D and d2 = e_1(α*) E(y). Both are lem:modrel, proved in Proposition 2.12 (the first is its first relation, not a definition). With c1 = C1cf = A_cf, and Cf taken **by convention** to be the identity on G_{T_F} (well defined by (C4); BS sec:cfgen describes Cf only through τ -> -1/(3τ), which does not fix the Hauptmodul), they are Cf(R_{ai}). A general solution of L2b gives the c1-polarisation, not Cf. BS's A satisfies L2b because A = θ of a period of L2b∘θ.

Theorem 2.7 identifies Definition 1.3 with the T_F limit of BS's canonical rows (R2:Definition 7); the z-only factor is exactly 1. Since round 5 that identification is **not** used in the proofs of Theorems LC_P2 and R_P2 or of Corollary SV_P2: everything in s.2.6 and s.3-s.6 is algebra on the closed form, including symplecticity (Proposition 2.10). The identification is needed for Proposition 2.11 and Proposition 13.2_P2.

**Pairing.** ⟨a, b⟩ := Sum η^{ab} a_a b_b, with η^00 = ε1 ε2 E^3/9, η^11 = -3 ε1 ε2 E, η^12 = η^21 = -3 ε1 ε2, all others 0. Note η^11 = E η^12 and η^22 = 0. In Theorems LC_P2 and R_P2 this is a definition. It agrees with lp2_lib etainvF and, for the last three values, with BS line 5395; matching it with eq:etaflat is needed only for the interface (s.2.7).

### 1.4 The factors
- **Vertex weight:** Delta_i.
- **Dilaton factor:** Dil^_i(d) := (-1)^d [z^{d-1}] P^_0^{(i)} for d >= 1, equivalently q_t^{(i)}(ψ) := Sum_{d>=1} Dil^_i(d) ψ^d = -ψ S(v_i, -ψ). Since [z^0]S = 1, we have Dil^_i(1) = -1 at every t: the d = 1 term is the dilaton shift, and the decorations of (0.1) have d >= 2. **Dil^_i uses only component 0 of the row.** So it does not depend on c1, and neither Lemma G nor the sign of c1 is involved. We write q_t := q_t^{(0)}.
- **Edge factor:** with edge variables ψ, w (rows at z = -ψ, z = -w; s.1.1) and P_ij(ψ, w) := -⟨P^^{(i)}(-ψ), P^^{(j)}(-w)⟩,

      Ehat_ij(ψ, w) := [P_ij(ψ, w) - P_ij(-w, w)]/(ψ + w),   Ehat_ij(m, n) := [ψ^m w^n] Ehat_ij.        (1.1)

  By symplecticity (Proposition 2.10), P_ij(-w, w) = -⟨P^^{(i)}(w), P^^{(j)}(-w)⟩ = -δ_ij Delta_i, which does not depend on w. So Ehat_ij(ψ, w) = [δ_ij Delta_i - ⟨P^^{(i)}(-ψ), P^^{(j)}(-w)⟩]/(ψ + w), which is the Givental edge in Delta-normalised form. Since η is symmetric, this gives the **edge symmetry** Ehat_ij(ψ, w) = Ehat_ji(w, ψ) for all i, j, so T(Gamma, c) does not depend on the orientation of the edges. (1.1) is a definition, and the conifold Lemma 1.1 (LSP) is the rank-free algebra behind it; see Lemma 1.4 for the code.

**Lemma 1.4 (code bookkeeping).** [P, by reading hodge_vertex/hv_lib3.m lines 57-60 (hodgeClass), 72-73 (hodgeNormU, usubOf), 75-78 (hCoeffs), 119-130 (vertexOptionsH), 131-138 (prepareGraphH), 142-143 (rowsTimes, evalPointRows); hodge_vertex/hv_lib4.m lines 49, 114 (hodgeNormU4, vertexOptionsH4, at genus 4); lp2_lib.m lines 37-64 (evalPoint, graphValue); lp2_explore/lp2_percolour.wls lines 41, 47, 50, 81-84, 105-106, 116-118.] The value computed at frame-B rows by graphValue, and by the separate series loop of lp2_percolour.wls lines 116-118, which has the same form, equals |Aut|^-1 Sum_c T(Gamma, c). (graphs.m vertexOptions, line 46, has pure ψ weights wk and is not the Hodge vertex.)

*Proof.*
- Vertex integral. hodgeClass[g] expands prod_a Sum_j (-1)^j λ_j u_a^j, and with u_a = 1/w_a, w = {E, -ε1, -ε2} (usubOf; lp2_percolour.wls lines 47, 106) this is Omega_w of s.1.5. hodgeNormU[g, a] is Int Omega_w ψ^a. vertexOptionsH produces the entries {hp, lp, dils, wgt} with wgt = hodgeNormU/(legs! nd!), dils the ordered tuples of degrees >= 2, and hp of total degree <= the budget (mixed degree). So at legs = 0 and nd = m_v, wgt is (1/m_v!) Int Omega_w ψ^a ψ^d as in (0.1).
- The code's vertex factor contains sqD^(2g_v - 2 + n_v + legs + #dil), with legs = 0 (G_g has no external markings; the code's gr["m"] is the number of legs, not the m_v of (0.1)) and sqD = Delta^{1/2}. This is Delta^{g_v-1}, times one sqD per half-edge, times one sqD per dilaton.
- One sqD per dilaton turns Dil[i][n] = (-1)^n Pnum[0,i][[n]] into Dil^_i(n).
- One sqD per half-edge turns Ed into the normalised edge. Indeed pre(m, n) = δ_ij δ_m0 δ_n0 - (-1)^{m+n} P_m η P_n has generating function N_ij(ψ, w) = δ_ij - ⟨P(-ψ), P(-w)⟩, and Ed(m, n) = Sum_{k=0}^{n} (-1)^k pre(m+k+1, n-k) has generating function [N_ij(ψ, w) - N_ij(-w, w)]/(ψ + w). The telescoping is ψ^M - (-w)^M = (ψ + w) Sum_{k<M} (-1)^k ψ^{M-1-k} w^k, and M = 0 drops out, so the Kronecker term never enters. Multiplying by sqD_i sqD_j gives (1.1).
- The frame-B twist is rowsTimes[hCoeffs[wB, -1, KMAX], KMAX] (lp2_percolour.wls line 50), i.e. multiplication by h_w^{-1}.
- The Leaf factors Lf occur only with markings, so not in G_g. ∎

This lemma relates (0.1) to the code. The GW meaning of Dil^ (Givental's T(z) = z(1 - R^{-1}(z) 1)) is **not** derived here; it belongs to Proposition 13.2_P2. Theorems LC_P2 and R_P2 are statements about the decoration sum as defined.

**Remark (signs).** [P] In the normalised form (0.1) no square root occurs. Delta_i = 3 ε1 ε2 E/v_i^2, and P^_a = Delta_i^{1/2} P_a^B does not depend on the sign of Delta_i^{1/2} because P_0^B carries Delta_i^{-1/2}, so that P^_a is built from S, θ, Y_i, C1i, d2 and θ log Delta_i, all single valued. Hence each T(Gamma, c) is well defined on any cover. The compatibility of branches in the code is built in, since sqDsym is transported with the rows. (The sign of c1 is a different matter: see Lemma G.)

### 1.5 Hodge side
For a vertex of genus h, the class is the **normalised** class

    Omega_w := prod_{a=1}^{3} Lambda^vee_h(w_a)/w_a^h,   w = (E, -ε1, -ε2),   Lambda^vee_h(u) := Sum_i (-1)^i c_i(E_Hodge) u^{h-i}.

(This is the convention of LSP Prop 13.2, of the code hv_lib3.m lines 2-3 and 57-60, and of the conifold's hvc_lib.m lines 3, 139. The scale law of s.7.4 depends on it: with the unnormalised product the exponent there would change.) It is a pull-back under forgetful maps, so the string and dilaton equations hold for h >= 2.

Definitions, for h >= 2. Let 𝐑 := C[[t_c]] (or Q(ε1, ε2)[[t_c]], or any t_c-adically complete Q-algebra of coefficients).

- **Base-point derivatives.** For g_1, ..., g_s in 𝐑[[ψ]], put

      𝔡^s(g_1, ..., g_s) := Sum_m (1/m!) Int_{Mbar_{h,s+m}} Omega_w prod_i g_i(ψ_i) prod_{k=1}^m T(ψ_{s+k}).

  This is a finite sum, since T = O(ψ^2) forces m <= 3h - 3 + s, and each power series is cut off by the dimension.
- **Fpot on the neighbourhood 𝔇 of q_0.** Let 𝔇 := q_0 + t_c 𝐑[[ψ]]. For q = q_0 + δ in 𝔇 and f_1, ..., f_n in 𝐑[[ψ]], define

      d^n Fpot_h(q)[f_1, ..., f_n] := Sum_{r>=0} (1/r!) 𝔡^{n+r}(f_1, ..., f_n, δ, ..., δ),   Fpot_h(q) := d^0 Fpot_h(q),

  which converges t_c-adically because δ = O(t_c). Note that δ may have a ψ-linear term: 𝔇 contains ρ q_t, whose ψ-coefficient is -ρ. So q need not be of the form -ψ + O(ψ^2).
- **The dilaton-shifted form.** If q = -ψ + O(ψ^2), put Fpot_h(q) := Sum_m (1/m!) Int_{Mbar_{h,m}} Omega_w prod_{i=1}^m (q + ψ)(ψ_i), a finite sum, and similarly for d^n Fpot_h(q)[f_1, ..., f_n]. On the overlap with 𝔇 the two definitions agree: expand prod (T + δ) binomially, and the factor (1/m!) C(m, r) reindexes to (1/r!)(1/(m-r)!). This form is what occurs at the regular colours, where q_t^{(j)} is not in 𝔇 and no rescaling or Taylor expansion is applied.
- For a multiset R = {j_1, ..., j_r} of integers: Xi_h(R) := d^r Fpot_h(q_0)[λ_{j_1}, ..., λ_{j_r}] = 𝔡^r(λ_{j_1}, ..., λ_{j_r}) = Sum_m (1/m!) Int_{Mbar_{h,r+m}} Omega_w prod_i λ_{j_i}(ψ_i) prod_k T(ψ_{r+k}).
- J_R := Sum j_i, D(R) := Sum (j_i - 1).
- Z[h, n] := Xi_h(1^n), so Z[h, 0] = Fpot_h(q_0).

**M_irr(h):** Xi_h(R) = 0 for every **non-empty** multiset R of integers **>= 2** with D(R) <= 2h - 3. The empty multiset is excluded, because Xi_h(∅) = Z[h,0] is the gap. By Theorem VI it is a non-zero element of Q(ε1, ε2); this does not specialise (s.7.7(f)), and no proof step uses it. Parts 0 and 1 are not in M_irr: part 1 is dealt with by the string identity, and part 0 never occurs (s.7.2).

The vertex factor. By the definition of Fpot_h, and because the decorations of (0.1) have d >= 2 with weight 1/m!, the sum over dilaton decorations at a vertex v of genus h >= 2, colour i, with half-edge series f_1, ..., f_n (f_r(ψ) = Sum_a (edge coefficient) ψ^a), equals

    Delta_i^{h-1} d^n Fpot_h(q_t^{(i)})[f_1, ..., f_n].       (1.2)

At colour 0, q_t is in 𝔇 by Lemma D0(iv). Edges do not factorise in general: a loop, or an edge between two rescaled vertices, is a tensor. (1.2) is multilinear in the half-edge slots and is applied term by term, leaf by leaf in the expansions of E00(i) and (M-b), and coefficient by coefficient in the edge variables.

### 1.6 Frames and modes
Frame B means R' = R_Z h_w^{-1}. Mode B means w = (E, -ε1, -ε2). The transfer of Lemma A holds for mode B only: there the leaf scale E produced by F2 equals the weight w_1, which is what allows Lambda^vee(w_1) to absorb the K-class (proof_status_mixed_scale.md s.1.1(c)). In mode A (w -> -w) the leaf scale is not a weight.


## 2. Rows on T_F

The closed form of Definition 1.3 is an algebraic object: all of s.2.1 (PF, Theorem A) and all of s.2.6, **including symplecticity** (Proposition 2.10), are statements about it alone. s.2.2-2.5 and s.2.7 identify it with BS's rows. They are not used in the proofs of Theorems LC_P2 and R_P2, and they are needed for Proposition 2.11 and Proposition 13.2_P2.

### 2.1 Picard-Fuchs for the closed form

Put a := 1 + E/z, e_n := V^n/(a)_n (n in Z, V = E v_i), with (a)_n the Pochhammer symbol. Then v S = z^{-1} Sum_{n>=1} r_{n-1} e_n. L_pf := prod_j (zθ + α_j) + z^3 y prod_{j=0}^{2} (3θ + j) (eq:pf, label at BS 2628, displayed about 2627-2632), and 𝒫 := e^{-u_i/z} L_pf e^{u_i/z} with θ u_i = Y_i (the operator of s.2.4 at T_F).

**Theorem A (R3:Theorem A).** [P] Let 𝗌 e_n = e_{n-1}, 𝖭 e_n = n e_n, 𝗍 = 3𝗌 - 𝖭. Then 𝒫 = (z^3/27)(𝗍^3 - 27 𝗌^3), and

    𝒫(Sum b_n e_n) = -(z^3/27) Sum_m [27(m+1) b_{m+2} - 3(3m^2+3m+1) b_{m+1} + m^3 b_m] e_m.

Hence Delta_i^{-1/2} h_w(z) S(v_i, z) e^{u_i/z} solves eq:pf identically, at i = 0, 1, 2.

*Proof.* We have, at every colour, θV = 1 - V/3, Y + α = E/V and y(V-3)^3 = 1.
- Conjugation gives 𝒫 = 𝖳^3 + 27y prod_{j=0}^{2} (𝖳 - E/3 + jz/3) with 𝖳 = zθ + Y + α.
- V e_n = (a+n) e_{n+1} and z(a+n-1) = E + nz give 𝖳 e_n = z e_{n-1} - (nz/3) e_n, i.e. 𝖳 = (z/3)𝗍, and 𝖳 - E/3 + jz/3 = -(𝖮 - jz)/3 with 𝖮 := E + z𝖭 - 3z𝗌.
- Let M be multiplication by z(V-3). Then 𝖮 = M𝗌 and [𝗌, M] = z, so M^2𝗌^2 = 𝖮(𝖮 - z) and M^3𝗌^3 = 𝖮(𝖮 - z)(𝖮 - 2z). Since yM^3 = z^3, 𝒫 = (z^3/27)(𝗍^3 - 27𝗌^3).
- 𝗍^3 e_n = 27e_{n-3} - 27(n-1)e_{n-2} + 3(3n^2-3n+1)e_{n-1} - n^3e_n gives the displayed coefficient. With b_n = r_{n-1} it vanishes by (Rec) for m >= 0, and trivially for m <= -1.
- Rigour (z-adic). Since 1/(a)_n = prod_{j<=n} z/(E + jz), each e_n with n >= 0 lies in V^n z^n Q[[z]], so Sum b_n e_n converges z-adically. 𝒫 has coefficients polynomial in z and never lowers the z-order. So [z^k] 𝒫(vS) is a finite expression, a Laurent polynomial in V and (V-3), and it vanishes in Q(t_i). This is the form used by Proposition 2.10.
- Rigour (V-adic, at the conifold). The same computation holds coefficientwise in V, with coefficients rational in z. On the V-adic side, e_n lies in V^n Q(ε)(z), 𝒫 e_n is the finite combination of e_{n-2}, e_{n-1}, e_n given by the displayed coefficient, and y = (V-3)^{-3} is analytic at V = 0. So 𝒫(vS) = 0 holds in Q(ε)(z)[[V]]. This is the form used by Z-infinity (Corollary 2.9).
- Finally, Delta_i^{-1/2} = v_i/sqrt(3 ε1 ε2 E), so Delta_i^{-1/2} S is a constant times v S = z^{-1} Sum b_n e_n, and the z-only factors h_w(z), z^{-1} commute with 𝒫. (R3 said "Delta_i^{-1/2} commutes with 𝒫"; it does not, as it depends on y. It is absorbed in v S, as here.) ∎

**Corollary 2.2 (R3:Cor A5, uniqueness).** [P] If b_n in C((z)), b_n = 0 for n <= 0 and 𝒫(Sum b_n e_n) = 0, then b_n = b_1(z) r_{n-1}, so Sum b_n e_n = z b_1(z) v S. Thus PF fixes the row only up to a factor depending on z alone. (Round-3 fix: the factor is z b_1(z), not b_1(z).)

### 2.2 Formal stationary phase (SP1-SP5)

Not used in the proofs of Theorems LC_P2, R_P2. The items were called F1-F5 in R3 and earlier rounds; they are renamed here because F1-F3 are also parts of Theorem LC_P2.

Part C of R3, in n dimensions. [S: written out, all formal.]

**The calculus.** Let A be the ring of sums Sum c_{k,β} z^k ξ^β with finitely many terms having k + |β|/2 <= N, for each N. Put ⟨f⟩ := e^Δ f|_{ξ=0}, with Δ = (z/2) H^{ab} ∂_a ∂_b. This is continuous on A. Define Λ(γ) := (det H)^{-1/2} ⟨e^{-φ^ng/z} γ⟩. Asym denotes this normalised formal Laplace expansion (BS line ~4245) at a non-degenerate critical point. The calculus is written for e^{-φ/z}; BS's e^{𝒲/z} is φ = -𝒲. In dimension one, SP5 with z -> -z is how the factors Θ(-z/a) arise; the difference of Gaussian conventions is the universal constant c of Proposition 2.6.

- **G1 (Gaussian integration by parts).** ⟨(Hξ)_a f⟩ = z⟨∂_a f⟩, because [Δ, ξ^b] = z H^{bc} ∂_c commutes with Δ.
- **SP2 (total derivatives vanish).** Λ(∂_a B - z^{-1}(∂_a φ)B) = 0, by G1.
- **FP (fixed critical point).** If the critical point sits at ξ = 0 **for every p** (the linear term of φ_p vanishes identically in p) and the critical value does not depend on p, then d/dp Λ_p(γ_p) = Λ_p(γ_p' - z^{-1} φ_p' γ_p). Proof: apply G1 twice to reproduce -(1/2) tr(H^{-1}H') and Δ'. Termwise legitimate by Lemma 2.3.
- **SP3 (moving critical point).** Translate the critical point to 0 and apply FP; the x_0' terms die by SP2. Result: d_p[e^{-φ_0/z} Λ] = e^{-φ_0/z} Λ(∂_p γ - z^{-1}(∂_p φ)γ). SP3 holds for any analytic amplitude γ.
- **SP1 (coordinate invariance).** With the branch convention (det H')^{1/2} = det Dg(0)·(det H)^{1/2}: the linear case is direct. In general, interpolate by g_s(ξ) = x_0 + s^{-1}(g(sξ) - x_0). Then H_s is constant in s, and X_s = (Dg_s)^{-1} ∂_s g_s vanishes at 0, with φ_s' = X(φ_s), γ_s' = X(γ_s) + div(X)γ_s. By FP and SP2, Φ'(s) = 0. Each z-coefficient is polynomial in s, so Φ(1) = Φ(0).
- **SP4 (Fubini).** If the phase splits exactly, Δ splits into commuting parts. The amplitude need not split. (This also disposes of round-1's "3d vs iterated Laplace" worry: in eq:asympRp2 the Hessian is block diagonal and the non-Gaussian part splits.)
- **SP5 (formal Stirling; formal version of BS exmpl:stirling, 4260, which BS state analytically for α in R^+).** For α != 0, Asym_{q=α}[e^{-(q - α log q)/z} q^{β-1} dq] = α^{β-1/2} Θ(z/α, β), with the branch of α^{β-1/2} fixed by the SP1 convention; this matters at α = -E. Scaling gives the form α^{β-1/2} F(w, β), w = z/α. SP2 with B = q^β gives (S): F(β+1) = (1 + βw)F(β). SP3 in α gives (D): ∂_β F = (β - 1/2) w F - w^2 ∂_w F. Writing log F = Sum g_k w^k, (S) makes the difference of two solutions constant in β, and (D) at index k+1 kills those constants, so the solution is unique. Θ satisfies (S) and (D), using B'_{k+1} = (k+1)B_k.

### 2.3 Critical points, analyticity, PF under the integral

Superpotential (eq:spotP2, 4301), in log coordinates x_j = log q_j: 𝒲 = Sum_j (e^{x_j} - α_j x_j) + y^{-1/3} e^{(x0+x1+x2)/3}, dΩ = (1/3) d^3x, P := y^{-1/3} e^{Sum x/3}. L^(Y) := prod_j (Y + α_j) + 27 y Y^3 (eq:LLMM, 2607).

**Lemma 2.1 (R2:Lemma 1).** [P]
- The critical points of 𝒲 are exactly the points with q_j = α_j + Y, P = -3Y, where L^(Y) = 0 and Y(α_j + Y) != 0.
- det Hess 𝒲 = W(Y)/3, and on L^ = 0: W = -Y L^'(Y).
- At T_F the three roots Y_i are finite, non-zero, simple and non-degenerate iff y is not in {0, -1/27}.
- At y = -1/27 only the colour-0 root escapes on the real branch, since 1 + 3ω^i t -> 1 - ω^i != 0 for i = 1, 2.
- θ u = Y at each critical point.

*Proof.* ∂_j 𝒲 = q_j - α_j + P/3. Put Y := -P/3; the constraint P^3 = prod q/y becomes L^(Y) = 0. Lifts to log coordinates differ by 2πiZ and have identical jets. Hess = diag(q) - (Y/3)J, J the all-ones matrix, and with e_2(q) = e_2(α) + 2e_1(α)Y + 3Y^2 this gives det = prod q - (Y/3)e_2(q) = W/3. From 81 y Y^2 = -3 prod(Y+α)/Y we get Y L^' = Y e_2(Y+α) - 3 prod(Y+α) = -W. On T_F, (Y+α)/Y = -3t_i. Finally θu = θ𝒲|_cr = -P/3 = Y. ∎

**Lemma 2.3 (R2:Lemma 2).** [P] On a parameter domain where the critical point x_c(p) is non-degenerate and analytic, and with an analytic branch of (det H)^{-1/2} on a simply connected set, every z-coefficient of Asym is analytic in p. For 𝒲, on a simply connected neighbourhood U of the semisimple locus in (y, α), the z^n coefficient of 𝒜 := Asym_{x^cr}[e^{𝒲/z} dΩ] lies in W^{-1/2} Q[α, Y][W^{-1}]. *Proof.* The z^n coefficient is the finite sum over m - k = n, 3k <= 2m, of (2^m m! k!)^{-1} (H^{ab}∂_a∂_b)^m [(-f_ng)^k g]|_{x_c}, polynomial in jets and in adj(H)/det H. The jets of 𝒲 of order |β| = m >= 2 are q_j δ_{β = m e_j} - Y/3^{m-1}. ∎

**Lemma 2.4 (R2:Lemma 3, with step 2 written out; R3:Cor B step (i)).** [P, using SP2, SP3] For every (y, α) in U, including α = α*, e^{u/z} 𝒜 is annihilated by L_pf.

*Proof (complete; this was the routine-but-unwritten item).* For an amplitude h(x, y) analytic near the critical point put I[h] := e^{u/z} Asym[h e^{𝒲/z} dΩ]. All amplitudes below are polynomials in q_j and P.

Two rules follow from SP3 and SP2. Since θ_y 𝒲 = -P/3 and ∂_j 𝒲 = q_j - α_j + P/3:

    (A)  zθ I[h] = I[zθ h - (P/3) h]                 (SP3 with p = log y, any analytic amplitude),
    (B)  I[z ∂_j h + (q_j - α_j + P/3) h] = 0        (SP2).

Adding (A) to (B) and to α_j I[h] gives, for each j,

    (C)  (zθ + α_j) I[h] = I[q_j h + z ∂_j h + zθ h].

Here θ acts on the explicit y-dependence of h at fixed x, and it commutes with ∂_x.

Step 1 (from Π(zθ + α_j) to q_0 q_1 q_2). Apply (C) three times, with amplitudes 1, q_2, q_1 q_2. Since ∂_1 q_2 = 0, ∂_0 (q_1 q_2) = 0 and θ q_j = 0 at fixed x:

    (zθ+α_2) I[1] = I[q_2],   (zθ+α_1) I[q_2] = I[q_1 q_2],   (zθ+α_0) I[q_1 q_2] = I[q_0 q_1 q_2].

So prod_j (zθ + α_j) I[1] = I[q_0 q_1 q_2]. This is SP3 applied three times and SP2 applied three times, with no division.

Step 2 (from q_0 q_1 q_2 to the third-order operator). θ P^m = -(m/3) P^m at fixed x. By (A) with h = P^m: zθ I[P^m] = I[-(mz/3) P^m - P^{m+1}/3], that is,

    I[P^{m+1}] = -z(3θ + m) I[P^m].

Hence I[P^3] = -z^3 (3θ+2)(3θ+1)(3θ) I[1], where the factors commute. This is SP3 applied three times.

Step 3. q_0 q_1 q_2 = y P^3, and multiplication by the scalar y passes through I by linearity. So

    prod_j (zθ + α_j) I[1] = y I[P^3] = -z^3 y prod_{j=0}^{2} (3θ + j) I[1],

which is L_pf I[1] = 0. ∎

Remark. In the form of round 2 (and of R3:Cor B) step 2 reads "iterating gives Π(zθ+α_j) ε = q_0 q_1 q_2 ε modulo total derivatives". The referee's concern was that the amplitudes in the intermediate steps are θ-derivatives and total derivatives; the formulation (C) avoids this, because the amplitude is explicit at each stage. The citation BS eq:oscintP (4315) is an alternative route at generic α (allowed), but it is not needed.

### 2.4 Uniqueness and normalisation transport

Let 𝒫 := e^{-u/z} L_pf e^{u/z} = Sum_{m=0}^{3} z^m 𝒫_m, with θu = Y.

**Lemma 2.5 (R2:Lemma 4).** [P] Let Ω' be a connected domain on which Y is holomorphic, L^(Y) = 0 and Y W(Y) != 0. Then:
- (i) 𝒫_0 = 0 and 𝒫_1 = -(W/Y)(θ + ½ θ log W). Hence ker 𝒫_1 = C·W^{-1/2}.
- (ii) If f = Sum f_k z^k, g = Sum g_k z^k have coefficients in O(Ω'), 𝒫f = 𝒫g = 0 and g_0 is not identically 0, then f = c(z) g for a unique c in C[[z]]. When α varies, c depends on α only, not on y.
- (iii) If moreover f, g extend analytically to a boundary point v0 of Ω', with g(v0) invertible in C[[z]], then c = f(v0)/g(v0). (Added hypothesis, round-2 referee: v0 = conifold is outside the good domain.)

*Proof.* (i): (Y + zθ)^n F = Y^n F + z[n Y^{n-1} θF + C(n,2) Y^{n-2} (θY) F] + O(z^2); use L^' = -W/Y and θY = -Y prod(Y+α)/W. (ii): the z^{k+1} coefficient is 𝒫_1 f_k + 𝒫_2 f_{k-1} + 𝒫_3 f_{k-2} = 0, and 𝒫 commutes with series in z alone; induct. (iii): evaluate. ∎

**Proposition 2.6 (R2:Prop 5).** [P mod (C1), (C2); uses SP1-SP5] Fix generic α, and let R^α_{0i} be BS's canonical row (eq:canR; (C1): e^{u_i/z} R^α_{0i} formally solves eq:pf; the citation is eq:asymp with eq:canR, not (Flatness)). Continue R^α_{0i} along any path in the semisimple locus (a germ along a path: R^α is defined near y = 0, and the monodromy around y* = u/(27(1-u)) permutes colours on the δ-family, so the statement is about germs, not single-valued holomorphy). Then R^α_{0i} = 𝒜(y, α, Y)/(c N(z)), where c in C^× is a universal constant with c^2 = ±1, fixed by the Gaussian square-root convention. In particular N depends neither on α nor on i.

*Proof.*
1. Near y = 0, 𝒜_i = κ_i(z) R^α_{0i} with κ_i independent of y, by (C1), Lemma 2.4 and Lemma 2.5.
2. Direct limit (this replaces the circular citation of BS 4337-4343, which BS deduce from eq:normNLR). At colour i, q_i = O(y). Put x_i = log y + log p, which adds the constant -α_i log y to u. Then 𝒲_y = yp - α_i log p + Sum_{j != i} (q_j - α_j x_j) + (p q_j q_k)^{1/3}, analytic at y = 0.
3. Set w' := (p q_j q_k)^{1/3}, a linear change in log coordinates (SP1), with dΩ = dlog w' dx_j dx_k. At y = 0, 𝒲_0 = Sum_{j != i} (q_j - (α_j - α_i) x_j) + (w' - 3α_i log w'), which decouples into three non-degenerate one-dimensional critical points (SP4).
4. By Lemma 2.3 at y = 0 and SP5: lim_{y->0} 𝒜_i = c prod_{a in {3α_i, α_j-α_i, α_k-α_i}} a^{-1/2} Θ(-z/a).
5. By (C2) (e^{-c_i} = prod_{k=0}^{4} Θ(-z/w_i^(k)) over the five tangent weights {3α_i, α_j-α_i, α_k-α_i, -ε1, -ε2}, Psi_{0i} = Delta_i^{-1/2} from Psi = Phi Delta^{-1/2}, BS 3014, eq:phi 3031), R^α_{0i}(0) = prod_5 w^{-1/2} Θ(-z/w).
6. The ratio, using Θ(-z/(-ε)) = 1/Θ(-z/ε), is κ_i = c sqrt(ε1 ε2) Θ(-z/ε1) Θ(-z/ε2) = c N.
7. Along any path both sides remain WKB solutions with the same Y, so Lemma 2.5 on a chain of discs keeps the ratio equal to cN. As N does not depend on i, relabelling by monodromy is harmless. ∎

So eq:normNLR (4327) is derived rather than cited.

**Theorem 2.7 (R2:Theorem 6, Prop 8, R3:Cor B; continuity and identification).** [P mod (C1)-(C3); uses SP1-SP5] Let y0 not be in {0, -1/27} and let Y* be a root at (y0, α*).
- (a) On a simply connected neighbourhood U of (y0, α*), every z-coefficient of R^α_{0i} extends analytically across the diagonal. Its value there is R^{T_F}_{0i} := [𝒜(·, α*)/(c N)]. Convergence is locally uniform with all y-derivatives. (N_i does not depend on α, so only the left side needs a limit.) This proves BS's bare "by continuity" assertion, about line 4345.
- (b) e^{u/z} R^{T_F}_{0i} solves eq:pf on T_F and its z^0 coefficient is Delta_i^{-1/2}, with the continued branch. The z^0 term fixes c.
- (c) The rows R_1, R_2 of eq:RHb, flatness, and Sum_{a,b} R_{ai}(z) η^{ab} R_{bj}(-z) = δ_ij pass to α = α*. This holds wherever A, D != 0, ∞, so not at y = -1/27. The coefficients of eq:RHb depend on α only through e_1(α), and θY = -Y prod(Y+α)/W is analytic on U.
- (d) For colours 1 and 2, (a), (b) also hold near y = -1/27, for row 0 (and for rows 1, 2 in the conifold frame only, since A has a logarithm there).
- (e) **Closed form.** R^{T_F}_{0i} = Delta_i^{-1/2} h_w(z) S(v_i, z). So P_0^{(i)} of Definition 1.3 is R^{T_F}_{0i}, with z-only factor exactly 1.

*Proof of (a)-(d).* Combine Proposition 2.6 with Lemma 2.3 and Lemma 2.4. The root stays finite and simple for (d).

*Proof of (e)* (BS 5216-5273, re-read; the reading of 5257-5267 was re-checked by a round-3 referee).
1. Colour 0, chart. The critical point is q^{1/3} = c'(1, ω, ω^2). This is a vector of q^{1/3}-coordinates (BS ~4370 writes it as if it were q; erratum, s.9); in q it is c'^3(1,1,1). Its entries are distinct, so s_n := e_n(q^{1/3}) is a local biholomorphism there (Disc = -27 s_3^2 != 0), and SP1 applies.
2. Colour 0, computation. In the s-chart 𝒲 = s_1^3 - 3 s_1 s_2 + (3 + y^{-1/3}) s_3 - E log s_3 splits exactly and the Hessian is block diagonal. The (s_1, s_2)-block yields Sum_l r^W_l (z/s_3)^l, which **defines** BS's Wick numbers r^W_l (BS 5252-5262). The s_3-integral is computed by the substitution s_3 = -q/V (q the variable of SP5), and then SP5 with α = -E, β = -l-1, where Θ(w, -m) = Θ(w)/prod_{j<=m}(1 - jw). This gives v^{l+1} Θ(-z/E)/prod_{j<=l+1}(1 + jz/E), up to a constant. **The sign (-1)^{l+1} from s_3 = -q/V cancels against (-E)^{-l-1}, so the prefactor carries no l-dependent sign.** The branch signs (at α = -E, and of square roots in SP1) go into a constant K_0, and r^W_0 = 1 is unaffected. Hence 𝒜_0 = K_0 Θ(-z/E) v S_Wick(v, z), with K_0 independent of z and y.
3. **The Wick numbers satisfy (Rec) (R3:Cor B).** By Lemma 2.4, e^{u/z} 𝒜_0 solves eq:pf. Θ(-z/E) is invertible and z-only, so 𝒫(v S_Wick) = 0, and Corollary 2.2 with r^W_0 = 1 gives r^W_l = r_l for all l.
4. Constants. N^{-1} Θ(-z/E) = h_w/sqrt(ε1 ε2), so R^{T_F}_{00} = K' h_w v S, with K' z-independent. The z^0 term must equal Delta_0^{-1/2} = v/sqrt(3 ε1 ε2 E) by (b), so K' = 1/sqrt(3 ε1 ε2 E). (No y-dependent factor can appear: it would not commute with eq:pf.)
5. Colours 1 and 2. t -> ωt amounts to x_0 -> x_0 - 2πi, which shifts 𝒲 by a constant and leaves dΩ unchanged. So the colour-i row is the colour-0 row with t replaced by ω^i t. ∎

Status of step 3: this makes the closed form S with r_l from (Rec) the same object as BS's Wick expansion. The dependence chain is SP1-SP5 (written) + Lemma 2.4 (written above) + (C1)-(C3).

Weights. Proposition 2.6 and Theorem 2.7 work over BS's generic weights. They are interface statements, and no specialisation of ε from them is used: Proposition 2.10 and everything in s.3-s.6 are proved directly for every complex ε with ε1 ε2 E != 0.

### 2.5 Status of the code [E only]
The stored file cPP_z9 is an implementation, checked against Definition 1.3: exactly at z^0..z^9 at 7 weight pairs; against eq:pf exactly through z^10 at 6 points; at z^0, z^1 symbolically by hand. The claim "genRmat computes the rows at all orders" is **not** proved: PF is solved at sample points; the rational reconstruction in u is heuristic; ansatz sufficiency (fixed denominators, bounded numerators) is unproved; and levels 0-6 come from a checkpoint (PFu6.m) whose run is not logged. pfcheck.pl cannot detect a wrong z-only factor by construction; c(z) = 1 for the stored rows rests on the round-1 row_vs_closed.pl comparison and the hand check of the z^0 and z^1 terms. No proof step uses the code.

### 2.6 Consequences of the closed form

**Corollary 2.8 (Lemma CL, R3 s.0, corrected; carrying identities to the conifold).** [P]
- At every colour, with c1 = C1cf, each z^n coefficient of P^_a^{(i)} is a **Laurent** polynomial in v_i (Lemma 4.1 has 1/v and 1/v^2) multiplied by z-free functions analytic at the conifold (C1i, θv, Δ_c, d2). The same holds in the large-radius frame (the factor h_w depends on z only), and for P_a^{B,(i)} = Delta_i^{-1/2} P^_a^{(i)} up to the constant sqrt(3 ε1 ε2 E)^{-1}. Note Delta_0^{±1/2} is proportional to v_0^{∓1}, so it is a factor of this Laurent type and is not analytic. So each such coefficient is a meromorphic germ at y = -1/27, and an identity between such series that holds on a punctured neighbourhood holds as an identity of Laurent expansions in t_c.
- **At colour 0 only, and in frame B only:** since v_0 = O(t_c), each [t_c^k] coefficient of ρ^e P^^{(0)} involves only finitely many v-orders. So [z^n] and [t_c^k] commute, and each [t_c^k] coefficient is a rational function of z. This fails at colours 1, 2, where v_j(0) != 0 and [t_c^k] of S(v_j, z) is only a formal power series in z, and it fails in the large-radius frame, since h_w is not rational. s.3.1 and Proposition H2 use only the colour-0, frame-B case.

**Corollary 2.9 (R2:Cor 9).** [P, from the closed form alone.] With the notation of s.1, σ_0 := ρ S(v_0, -ψ):
- (F2) S(0, z) = E/(E + z): only l = 0 survives at v = 0.
- (L1) Write S(v_0, z) = Sum_l r_l (vz)^l / prod. The l-th summand of P_0^{B,(0)} = v S/sqrt(3 ε1 ε2 E) carries v^{l+1} and the proper fraction z^l/prod_{j<=l+1}(1 + jz/E). Its partial fractions lie in span{λ_1, ..., λ_{l+1}}. As v = t_c·(z-free unit):
  - [t_c^k] P_0^{B,(0)} lies in span{λ_1, ..., λ_k};
  - [t_c^k](ρ S) lies in span{λ_1, ..., λ_{k+1}} (one more, because of the division by t_c).

  These are two different objects, and the round-2 referee asked for them to be distinguished. Neither has a constant or a polynomial part.
- (Form) [z^k](v h_w S) is a polynomial in v of degree <= k+1, equivalently of the stored shape N_k(Y)/(E+3Y)^{k+1}, deg N_k <= k. (Note this is not the statement [z^k] S has degree <= k in v, which is used for F3 below and is immediate from the definition of S.)
- (R) For i = 1, 2, every z-coefficient of P_0^{(i)}, P_1^{(i)}, P_2^{(i)} with c1 = C1cf is analytic at y = -1/27, in either frame. They are rational in t and regular where t(1 + 3ω^i t) != 0 and 1 - ω^i != 0; C1cf and XX are analytic with C1cf(0) = -1, and Δ_c d2 = E(C1i^2 - 1)/2.
- (Z-infinity) lim_{z->∞} z v S(v, z) = -t_c coefficientwise in v. Equivalently lim z ρ S = E.

*Proof of (Z-infinity).* Coefficientwise in v, f := v S = S_∞/z + O(z^{-2}), with S_∞ = Sum r_l (Ev)^{l+1}/(l+1)!. Write the conjugated operator as 𝒪 = zO_1 + z^2 O_2 + z^3 O_3, with top part O_3 = Δ_c θ^3 + 27 y θ^2 + 6 y θ = L2b∘θ. The V-adic form of Theorem A gives 𝒪f = 0 coefficientwise in v, as identities of rational functions in z. Taking z^{-2} 𝒪 f as z -> ∞ gives L2b(θ S_∞) = 0. The exponents at Δ_c = 0 are (0,0), so θ S_∞ = κ C1cf = κ θ t_c, and S_∞ = κ t_c because S_∞(0) = 0. The v^1 coefficients (E v against -E v, Lemma 1.2(b)) give κ = -1. ∎

(Z-infinity is exactly the place where the normalisation of t_c enters: it identifies S_∞(v_0) = -t_c with t_c normalised by C1cf = θ t_c and C1cf(0) = -1.)

**Proposition 2.10 (symplecticity, K = 1, orthogonality; algebraic proof, round 5).** [P] Let ε1, ε2 be complex numbers with ε1 ε2 E != 0, or indeterminates. For any non-vanishing section c1 with XX = θ log c1, and all i, j,

    Sum_{a,b} P_a^{(i)}(z) η^{ab} P_b^{(j)}(-z) = δ_ij.

The same holds in frame B since h_w(z) h_w(-z) = 1, so ⟨P^^{(i)}(z), P^^{(j)}(-z)⟩ = δ_ij Delta_i. By the drop-out S0 below, the left side depends only on P_0, and each of its z-coefficients is rational in t. So the identity holds at each z-order as an identity of rational functions of t, and (K=1) and (M-d) at the conifold follow by expanding in t_c, with no continuation or branch argument.

The proof uses only Theorem A and algebra on the closed form. It uses neither s.2.2-2.5, nor eq:RHb, (Symplecticity), lem:modrel or (C1)-(C3); and no genericity of ε, so no specialisation argument is needed.

*Proof.* Work over K := C(t), y = t^3, θ = (t/3) d/dt, with coefficients series in z. Write Th_i := Y_i + zθ on K[[z]], and at colour i put

    a := v_i S(v_i, z),   b := Th_i a,   c := Δ_c Th_i b.

Let 𝔅(u, w) := -(E^3/27) u_0 w_0 + E u_1 w_1 + u_1 w_2 + u_2 w_1, so that η^{ab} = -3 ε1 ε2 𝔅. Put Π̃_ij(z) := 𝔅((a, b, c)^{(i)}(z), (a, b, c)^{(j)}(-z)).

S0 (drop-out, for any c1). Put Q_1 := Th_i P_0^B. Using θ C1i = -C1i XX and Δ_c d2 = E(C1i^2 - 1)/2, Definition 1.3 gives P_1^B = C1i Q_1 and

    P_2^B = C1i^{-1} [𝒯 - z Δ_c XX Q_1] - C1i (E/2) Q_1,   𝒯 := Δ_c Th_i Q_1 + (E/2) Q_1.

In the pairing the C1i^2 terms from η^11 and η^12 cancel (using η^11 = E η^12), and so do the XX terms. **This requires the same c1 at both colours**, which holds since c1 depends on y only. There remains η^00 P_0^{(i)} P_0^{(j)'} - 3 ε1 ε2 (Q_1^{(i)} 𝒯^{(j)'} + 𝒯^{(i)} Q_1^{(j)'}), with ' meaning z -> -z. With c_0 := 3 ε1 ε2 E we have P_0^B = a/sqrt(c_0), Q_1 = b/sqrt(c_0), 𝒯 = (c + (E/2) b)/sqrt(c_0), and the h_w factors contribute h_w(z) h_w(-z) = 1. So the frame-B pairing equals -3 ε1 ε2 Π̃_ij/c_0 = -Π̃_ij/E, for every c1. It remains to prove Π̃_ij = -E δ_ij.

S1 (first-order system). Conjugating eq:pf (T_F, α = E/3) gives the stripped form Δ_c Th^3 + (E + 27yz) Th^2 + (E^2/3 + 6yz^2) Th + E^3/27, which annihilates a by Theorem A. Since Th(Δ_c f) = Δ_c Th f + 27yz f, the 27yz terms cancel, and

    Th a = b,   Th b = c/Δ_c,   Th c = -(E/Δ_c) c - (E^2/3 + 6y z^2) b - (E^3/27) a.

Write this as Th (a, b, c) = 𝖣 (a, b, c). The matrix 𝖣 depends only on y and z^2, not on the colour and not on the sign of z.

S2 (self-adjointness). For any vectors u, w,

    𝔅(𝖣u, w) = -(E^3/27)(u_1 w_0 + u_0 w_1) + u_2 w_2/Δ_c - (E^2/3 + 6y z^2) u_1 w_1,

because the two terms E u_2 w_1/Δ_c cancel. This is symmetric in u, w, so 𝔅(𝖣u, w) = 𝔅(u, 𝖣w).

S3 (ODE). Let u := (a, b, c)^{(i)}(z), w := (a, b, c)^{(j)}(-z). By S1 at z and at -z: zθ u = 𝖣u - Y_i u and zθ w = Y_j w - 𝖣w. Since 𝔅 has constant coefficients, S2 gives

    (Y_i - Y_j + zθ) Π̃_ij = 0.

S4 (i != j). Y_i - Y_j is a non-zero element of K and Π̃_ij lies in K[[z]], so induction on the z-order gives Π̃_ij = 0.

S5 (i = j). θ Π̃_ii = 0, so each z-coefficient is constant in t. Colour i is the image of colour 0 under the field automorphism t -> ω^i t, which fixes y, Δ_c and θ. Hence Π̃_ii = Π̃_00 =: κ(z) for all i.

S6 (the constant). Each z-coefficient of a lies in v C[v]. Next, b = (1 - Ev/3) S + z(1/E - v/3) d_v(vS), using Y = -E/3 + 1/v and θv = 1/E - v/3, so each z-coefficient of b lies in C[v]. And c = (Δ_c Y) b + Δ_c zθ b, with Δ_c Y = -(E/3) Δ_c + E t (1 - 3t + 9t^2). All are regular at t = -1/3, i.e. v_0 = 0. There a = 0; b = S(0, z)(1 + z/E) = 1, since S(0, z) = E/(E + z); Δ_c = 0; Δ_c Y = -E; so c = -E. Thus κ(z) = E·1·1 + 1·(-E) + (-E)·1 = -E, exactly in z. ∎

*Remarks.*
- Second route, through BS. At c1 = A, Theorem 2.7(c) and BS (Symplecticity) (3496) give the same identity on a punctured domain, and S0 extends it to every c1. This needs Definition 1.3 at c1 = A to be eq:RHb, i.e. Proposition 2.12, and also Theorem 2.7(e) to identify row 0 with BS's R_0. It was the proof in round 4 and is no longer used.
- Check, not a proof step: at v -> 0 on colour 0, P_0 = O(v), Q_1(0) = h_w/sqrt(3 ε1 ε2 E), 𝒯(0) = -(E/2) Q_1(0), and so the pairing is h_w(z) h_w(-z) = 1. Round 1's Lagrange identity and (K-ii) are not needed.
- This discharges the item marked [O, optional] in round 4 (an algebraic proof of Π_ij = δ_ij from the closed form). It has had one referee pass (sound with minor fixes, applied here; the referee re-derived Theorem A and S0-S6 by hand).

**Proposition 2.11 (R2:Prop 11; Givental-Teleman on T_F).** [P mod (C1)-(C3), (C5)] Not used in Theorems LC_P2, R_P2. For g >= 2, GW_g^{T_F} equals the graph sum built from R^{T_F} at every y not in {0, -1/27}, by continuation. *Proof.* Along eq:adelta (3232), GW_g(y, δ) = Sum φ_{k,m}(y, δ) F^k A^m (C5), with φ_{k,m} symmetric rational functions of the roots, φ_{k,m} = P_{k,m}/Q with coefficients in C{δ}[y]. If δ divided Q, some φ_{k,m}(y0, ·) would blow up, contradicting Theorem 2.7 (each graph term is analytic at (y0, α*)). Both sides of Q·GW_g = Sum P_{k,m} F^k A^m are regular at δ = 0 (the left by BS line ~4056), so the identity holds at δ = 0; the pointwise limit of the graph sum is the T_F graph sum. ∎ Note: BS line ~4056 is argued only briefly; see s.8.3. BS's prop:fingenR (3611) is stated for rational tori only; at T_F, R_0's independence of A follows from Theorem 2.7(e), since the row is rational in t.

### 2.7 lem:modrel and the interface with BS's rows

This subsection is not used in the proofs of Theorems LC_P2, R_P2 or Corollary SV_P2. It closes round 1's item N4 and is needed for the interface statement "Definition 1.3 at c1 = A is eq:RHb" in Proposition 2.11 and Proposition 13.2_P2.

**Proposition 2.12 (lem:modrel, BS 2695; both relations).** [P, no citation, no GW input, no finite-order step] Define A, D, E(y), X by eq:Ifun (2621) and eq:ADEFXdef (2681). Here E(y) is BS's function, always written with its argument to distinguish it from the weight E. Then, as elements of Q[[y]],

    A^2 D = X,     2 E(y) = D - X,

and eq:DEXrel (2702) and eq:limitconv (2711) hold. So these are identities of meromorphic germs on |y| < 1/27 (holomorphic where A != 0; no use is made of the location of zeros of A) and along any continuation.

*Notation.* a_d := (3d)!/(d!)^3, H_n harmonic numbers, Δ_c := 1 + 27y, X := 1/Δ_c, L2b := θ^2 + 3y(3θ+1)(3θ+2) = Δ_c θ^2 + 27yθ + 6y.

*Proof.*

Step 0 (read off eq:Ifun). With s = 1/z, the d-th summand (d >= 1) of I/(z y^{H/z}) is

    y^d prod_{m=0}^{3d-1} (-m - 3Hs) / prod_j prod_{m=1}^d (m + (H + α_j)s) = 3H (-y)^d ((3d-1)!/(d!)^3) s [1 + s(3H H_{3d-1} - (3H + e_1(α)) H_d) + O(s^2)].

Multiplying by z y^{H/z} and using 3d (3d-1)! = (3d)! gives

    A = Sum a_d (-y)^d,   G := θ I_1^{[1]}/e_1(α) = -Sum a_d H_d (-y)^d,   P_I := θ I_1^{[2]} = A log y + h,   h := Sum 3 a_d (H_{3d} - H_d) (-y)^d.

(For h: the term from I_0^{[1]} contributes (a_d/d)(-y)^d, and 1/d + 3H_{3d-1} = 3H_{3d}.) So E(y) = θ(G/A), D = θ(P_I/A), and E(y) does not depend on α. No reduction modulo the cubic relation occurs at this order. This agrees with eq:Ian.

Step 1 (L2b on series). For φ = Sum φ_d (-y)^d, [(-y)^d] L2b φ = d^2 φ_d - 3(3d-1)(3d-2) φ_{d-1}. Since d^2 a_d = 3(3d-1)(3d-2) a_{d-1}:
- (i) L2b A = 0;
- (ii) if φ_d = a_d g_d, then [(-y)^d] L2b φ = d^2 a_d (g_d - g_{d-1});
- (iii) L2b is injective on y Q[[y]], because the lowest term d >= 1 has coefficient d^2 φ_d != 0.

Product rule: L2b(A l) = Δ_c (2 θA θl + A θ^2 l) + 27y A θl.

Step 2 (key identity). Sum_d a_d (3H_{3d} - H_d) (-y)^d = -A log(1 + 27y).   (K)

Both sides lie in y Q[[y]]. For the left side, g_d = 3H_{3d} - H_d gives g_d - g_{d-1} = 3/(3d-1) + 3/(3d-2) = 9(2d-1)/((3d-1)(3d-2)), so by (ii) the (-y)^d coefficient of L2b is 27(2d-1) a_{d-1}. For the right side, θ log Δ_c = 27yX and θ^2 log Δ_c = 27yX^2, so the product rule gives L2b(A log Δ_c) = 54y θA + 27y A, whose (-y)^d coefficient is -27(2d-1) a_{d-1}. Hence L2b of the right side has the same coefficient. Apply (iii). Comparing coefficients in (K), after dividing by 27^d, is exactly eq:limitconv.

Step 3 (eq:DEXrel). 2G - h = -Sum a_d (3H_{3d} - H_d)(-y)^d = A log Δ_c by (K). So 2G - P_I + A log(y/Δ_c) = 0, which is eq:DEXrel.

Step 4 (second relation). Divide by A: G/A = (P_I/A - l)/2 with l = log(y/Δ_c). Since θ l = 1 - 27yX = X, applying θ gives 2E(y) = D - X.

Step 5 (first relation; this replaces the Zagier-Zinger citation). L2b P_I = 0: by (ii), with g_d = 3(H_{3d} - H_d), [(-y)^d] L2b h = 27(2d-1) a_{d-1} - 2d a_d. By the product rule, L2b(A log y) = 2Δ_c θA + 27yA, with coefficient 2d a_d - 27(2d-1) a_{d-1}. The sum is 0. Put W_r := A θP_I - P_I θA = A^2 D. Then θW_r = -27yX W_r, so θ(Δ_c W_r) = 0. W_r = A^2 + A θh - h θA is a log-free series with W_r(0) = 1. So Δ_c W_r = 1. ∎

**Dictionary.** [P] With d2 := E(C1i^2 - 1)/(2Δ_c) (s.1.2), the weight E, and c1 = A: Δ_c c1^2 = A^2/X = 1/D by the first relation, and d2 = (E/2)(X/A^2 - X) = (E/2)(D - X) = e_1(α*) E(y) by both relations. So Definition 1.3 at c1 = A has the recursion coefficients of eq:RHb. (Round 4 attributed c2 = D to eq:ADEFXdef. That was wrong: it is the first relation of lem:modrel.) The identification of row 0 is a separate statement, Theorem 2.7(e).

**Remark 2.13 (audit of BS; not used by Propositions 2.10, 2.12).** By reading BS_v2 main.tex; line numbers as of that file.
- (a) eq:RHb does not depend on lem:modrel. eq:SJ2 (2590) is pure GW axioms, eq:ssH1 (3077) is J = I, and the Divisor equation (about 3098) gives C^a_{1b} = A^{-1} θ s^a_b. So C^1_11 = e_1 E(y)/A and C^2_11 = D/A by eq:ADEFXdef alone, and C^0_11 = η^00 η_11 = 0 (commented pre-modrel forms survive at about 3110-3118). eq:SHb (display about 3299, label 3305) cites eq:allc, whose printed row 2 has modrel applied, but it uses only these unsimplified entries; eq:RHb follows from eq:canR. Suggested fix in BS: cite eq:ssH1 at 3299.
- (b) (Symplecticity) does not depend on lem:modrel. eq:sympl (3084) follows by differentiation in t, using only η-self-adjointness of the quantum product and S = 1 + O(1/z). BS 3543 passes to R in one clause; the formal passage is (Y_i - Y_j + zθ) B_ij = 0 for B_ij := Sum R_ai(z) η^ab R_bj(-z), which forces B_ij = 0 for i != j and θ B_ii = 0, and (Normalisation) with c_j odd in z gives B_ii = 1.
- (c) (In this remark e_k := e_k(α).) Where lem:modrel does enter BS downstream: printed eq:allc rows 2-4; eq:phi (Phi_2i through C^1_12, C^2_12); (Phi^{-1})_{i1}; eq:XFABC (2857); eq:Sain (the proof of prop:fingenS at about 3349) and prop:fingenR (3611); lines 3856 and 4836. BS's written derivation of Delta_i^{-1} = 1/(3 ε1 ε2 W_i) passes through these. A modrel-free derivation: at generic α, Phi_0i = 1, Phi_1i = Y_i/A, and idempotency with the values of (a) gives Phi_2i = Y_i(Y_i - e_1 E(y))/(AD) (D is not identically 0, as D(0) = 1). Then η_ab = Sum_i Phi_ai Phi_bi/Delta_i with a = 0, b = 0, 1, 2 gives Sum_i Delta_i^{-1} Y_i^k = δ_{k0}/(3 ε1 ε2 e_3) for k = 0, 1, 2, and Lagrange interpolation at 0 (the Y_i are distinct at generic α) with prod Y = -e_3/Δ_c gives Delta_i^{-1} = 1/(3 ε1 ε2 W_i), using W = -Y L^'(Y).
- (d) A GW proof of both relations, independent of Proposition 2.12, given the eigenvalues of eq:uili (2956). Y^3 = -(e_1 Y^2 + e_2 Y + e_3)/Δ_c gives Sum Y^3/Delta_i = -1/(3 ε1 ε2 Δ_c) and Sum Y^4/Delta_i = e_1/(3 ε1 ε2 Δ_c^2). Then η_12 = -1/(3 ε1 ε2) gives A^2 D = X, and η_22 = e_1/(3 ε1 ε2) gives e_1 X (X + 2E(y))/(A^2 D^2) = e_1, i.e. 2E(y) = D - X, for e_1 != 0 and D not identically 0. Since E(y) is α-free, one generic α suffices.
- (e) Round-4 critic 1's proposal at z^0 is correct. R_2^{RHb} = P_2|_{c1=A} + (1/D)(d2 - e_1 E(y)) R_1 (this uses A^2 D = X for the prefactor). The closed form at z^0 has pairing 1, from ε1 ε2 E^3/9 - 6 ε1 ε2 Δ_c Y^3 - 3 ε1 ε2 E Y^2 = 3 ε1 ε2 E (Y + E/3)^2 and Δ_c Y^3 = Y^3 - (Y + E/3)^3. So (Symplecticity) forces d2 = e_1 E(y), since R_1^{z^0} = Y Delta^{-1/2}/A != 0. The inputs are (a), (b), the modrel-free Delta_i of (c) or Lemma 2.5 with (Normalisation), and the first relation. It is superseded by Proposition 2.12.
- (f) Conceptual meaning: 2E(y) - D = -X is the trace identity for H*: C^1_11 + C^2_12 = e_1 (2E(y) - D)/A must equal Sum Y_i/A = -e_1 X/A. Likewise tr A_0 = -C1i E/Δ_c in Proposition 3.1, and the (1,2) symmetry entry of A_0 is exactly 2Δ_c d2 = E(C1i^2 - 1).

**Scope.** Propositions 2.10 and 2.12 reprove BS's (Symplecticity) on T_F independently of eq:sympl and of Theorem 2.7(c), but not independently of BS altogether: identifying row 0 of Definition 1.3 with BS's R_0 still needs Theorem 2.7(e), hence SP1-SP5, (C1) and (C2). The proof of Theorem 2.7(e) uses neither eq:sympl nor Theorem 2.7(c), so there is no circularity. BS lem:modrelorb (about 4581) is also proved only by a sketch; it is outside this chain and was not examined.


## 3. The singular colour: F1, F2, L1, Z-infinity, closure, K = 1, Dil_0

Throughout s.3-s.5, v := v_0, on the real branch, c1 = C1cf. Then θv = 1/E - v/3, Y = -E/3 + 1/v, θ log v = Y/E, Delta_0 ∝ v^{-2}.

### 3.1 Row 0 at colour 0
Put σ := ρ P^^{(0)}(z = -ψ) and σ^(k) := [t_c^k] σ. By Corollary 2.8 (the colour-0, frame-B case), the σ^(k) lie in the space V := C[ψ] ⊕ span{λ_j : j >= 1}; partial fractions make the decomposition unique. Write σ_0^(k) = Sum_j β_j^(k) λ_j. Translating Corollary 2.9 and Proposition 2.10:

- **(F1 at colour 0)** Lemma 1.2(b).
- **(L1 for σ_0)** σ_0^(-1) = 0 and σ_0^(k) lies in Λ_{k+1}.
- **(F2)** σ_0^(0) = E λ_1. [S(0, -ψ) = E/(E - ψ), ρ(0) = 1.]
- **(Z-infinity)** Sum_j β_j^(k)/j = E δ_{k0}. [lim_{z->∞} z ρ S = E, ψ = -z, and ψ λ_j -> -1/j.]
- **(K=1)** With ' meaning ψ -> -ψ: η^00 σ_0 σ_0' + η^11 σ_1 σ_1' + η^12 (σ_1 σ_2' + σ_2 σ_1') = 3 ε1 ε2 E^3/t_c^2 as an identity of Laurent expansions in t_c with coefficients in V ⊗ V. This is Proposition 2.10 at i = j = 0 and Lemma 1.2(b). Since Proposition 2.10 is an identity of rational functions of t at each z-order, no continuation or branch argument is needed; Corollary 2.8 gives the coefficients in V.

### 3.2 Dil_0: L1 and F3 restated (Lemma D0)

**Lemma D0.** [P, from Corollary 2.9]
- (i) ρ q_t = -ψ σ_0.
- (ii) [t_c^0](ρ q_t) = q_0 = E - E^2 λ_1.
- (iii) For k >= 1, [t_c^k](ρ q_t) = -E Sum_j (β_j^(k)/j) λ_j lies in Λ_{k+1}, with no constant.
- (iv) Each Dil^_0(d) is analytic at t_c = 0, with Dil^_0(d; 0) = -E^{1-d}.

*Proof.*
(i) By s.1.4, q_t = -ψ S(v, -ψ), and σ_0 = ρ S(v, -ψ).
(ii)-(iii) By (LAM1), [t_c^k](ρ q_t) = -ψ Sum_j β_j λ_j = Sum_j β_j/j - E Sum_j (β_j/j) λ_j. By Z-infinity, the constant is E δ_{k0}. At k = 0, β_1^(0) = E (F2), so the value is E - E^2 λ_1 = q_0.
(iv) Dil^_0(d) = (-1)^d [z^{d-1}] S. By the definition of S (s.1.3), [z^k] S is a polynomial in v of degree <= k, and v is analytic. At v = 0: (-1)^d [z^{d-1}] E/(E+z) = (-1)^d (-1/E)^{d-1} = -E^{1-d}. ∎

The vanishing of the constant for k >= 1 is exactly Z-infinity, i.e. the normalisation S_∞(v) = -t_c. A round-4 referee re-checked the period statement L2b(θ S_∞) = 0 exactly through V^14 [E; zinf.py]; the proof is in Corollary 2.9.

(Citation corrected from round 4: F3 for Dil_0 uses the definition of S, not (Form), which concerns [z^k](v h_w S).)

### 3.3 Closure (conifold-frame QDE). For general c1, not used in the proof of Theorems LC_P2, R_P2.

The c1 = 1 form of the closure is the first-order system S1 of Proposition 2.10, which is derived there directly from Theorem A. So Proposition 3.1 itself is not an input.

**Proposition 3.1 (round 1 closure; independently re-derived by the round-1 verifier and, for general c1 including the z^1 cancellation and self-adjointness, by the round-5 modrel author).** [P, from Theorem A] Let c1 be any invertible section, XX = θ log c1, and 𝐫 = (P_0^B, P_1^B, P_2^B) the vector of unnormalised rows of Definition 1.3 at any colour (the z-only factor h_w does not matter). With Th := Y + zθ,

    C1i Th 𝐫 = (A_0(y, C1i) - z^2 c1 Rd E_21) 𝐫,

where E_21 is the elementary matrix with entry (row 2, column 1), counting rows and columns from 0, and where
- A_0 = [[0, 1, 0], [0, C1i d2, C1i^3/Δ_c], [-E^3/27, -c1(E d2 + Δ_c d2^2 + E^2/3), -C1i(E/Δ_c + d2)]] is z-free, η-self-adjoint and independent of XX;
- Rd := Δ_c θXX + Δ_c XX^2 + 27y XX + 6y is the Riccati defect. c1 solves L2b iff Rd = 0.

So at c1 = C1cf the conifold-frame QDE z d/dt_c 𝐫 = (A_0 - C1i Y_i) 𝐫 holds at colour i, with A_0 z-free and self-adjoint.

*Proof.* Substitute P_0, Th P_0 = c1 P_1, Th^2 P_0 = P_2/(Δ_c c1) + (z XX + d2) c1 P_1, and reduce Th^3 P_0 with the stripped form of eq:pf, Δ_c Th^3 + (E + 27yz) Th^2 + (E^2/3 + 6yz^2) Th + E^3/27 = 0 (Theorem A). Collect powers of z: the P_2-coefficient is z-free since 27yz/Δ_c cancels. In the P_1-coefficient the z^1 term cancels (the identity θ d2 = -E XX C1i^2/Δ_c - (27y/Δ_c) d2 is one way to see it), and the z^2 term is -c1 Rd. The P_0-coefficient in row 2 is -E^3/27. Self-adjointness of A_0 uses η^11 = E η^12 and η^00 = -(E^3/27) η^12. ∎

**On-shell.** XX on the conifold section is the unique Riccati solution analytic at Δ_c = 0, XX(0) = -2/9, because the recursion (n+1) x_n = n x_{n-1} + [XX^2]_{n-1} + (2/9) δ_{n1} is non-resonant. Off shell, the z-even self-adjoint term E_21 drops out of the pairing (which is why symplecticity holds off shell, Proposition 2.10), but it produces the pole 3 ε1 ε2 E^2 (XX(0) + 2/9)/t_c of Ehat_00 [E]. (Round 1's ε1 ε2 E^2 (XX(0) + 2/9) is per 1/s, with t_c ~ 3s; it is the same statement.) The proof below does not use closure for general c1, or the flow. On-shell enters only through c1 = θ t_c, i.e. (I1), (I2) of Lemma 4.2, combined with Z-infinity; it is measured by H2(d).


## 4. Proposition H2

Everything in this section is all-orders algebra on the closed form, using only the inputs of s.3.1. It uses neither the closure lemma, the rank-3 flow nor round 1's (I1)/(I2) as hypotheses; these identities are proved here.

**Lemma 4.1 (closed forms for rows 1 and 2, any section c1; R3:Lemma 1).** [P]

    P^_1 = E C1i (θv/v) S^[1],
    P^_2 = E c1 Δ_c θv [E θv S^[2]/v^2 + (z θ log(C1i θv) - d2) S^[1]/v].

*Proof.* Use the identity Y + z θ log v = (E + z) θv/v, together with the termwise identities (E + (l+1)z)/prod_{j<=l+1} = E/prod_{j<=l} and (E + lz)/prod_{j<=l} = E/prod_{j<=l-1}, applied to the rows of Definition 1.3. ∎

**Lemma 4.2 (conjugation by ρ; R3:Lemma 2).** [P] Put g := θ t_c and D_g := Y + z g (1/t_c + d/dt_c). Then:
- (G) For any section: σ_1 = C1i D_g σ_0 and σ_2 = Δ_c c1^2 (D_g - d2) σ_1. This is conjugation by ρ, using θ log ρ = θ log v - θ log t_c and θ = g d/dt_c.
- On shell (c1 = g), with D_c := Ŷ + z/t_c + z d/dt_c:
  - (I1) σ_1 = D_c σ_0, where Ŷ := C1i Y = E d log v/dt_c = E/t_c + E/18 + O(t_c);
  - (I2) σ_2 = κ D_c σ_1 + μ σ_1, where κ := c1^3 Δ_c = -t_c + O(t_c^2) and μ := (E/2)(c1^2 - 1) = O(t_c).

Remark (not an input). This agrees with round 1's L_1 (here Ŷ, renamed to avoid a clash with the leaf statement L1), κ, μ.

**Lemma 4.3 (R3:Lemma 3).** [P] For F in t_c^{-1} V[[t_c]]:

    [t_c^n] D_c F = (E - (n+2)ψ) F^(n+1) + Sum_{m>=0} Ŷ^(m) F^(n-m),

and (LAM2) holds. Here Ŷ^(m) := [t_c^m] (Ŷ - E/t_c).

**Proposition H2.** [P; the inputs of s.3.1 are all P]
- (a) σ^(k) = 0 for k <= -2, and σ^(-1) = (0, E, -E^2).
- (b) For k >= 0, σ_0^(k) and σ_1^(k) lie in Λ_{k+1}, and σ_2^(k) lies in Λ_k. No component has a ψ-polynomial part or a constant at any k >= 0.
- (c) σ^(0) = (E λ_1, (E^2/9) λ_1, 0).
- (d) Within the family c1 = C1cf f(Δ_c) (f analytic), the ψ-linear coefficient of σ_2 is E Δ_c C1cf^2 θ(f)/t_c. It vanishes iff f is constant. Analyticity of d2 then forces f = ±1. So in this family, the absence of ψ-polynomial parts selects the conifold section, up to sign. (The statement is confined to this family.)

*Proof.*

**Step 1 (σ_1).** Apply (I1) and Lemma 4.3 with F = σ_0.
- n = -2 and below: the coefficient is 0 since σ_0^(k) = 0 for k <= -1.
- n = -1: (E - ψ) E λ_1 = E, by (LAM2) and F2. So σ_1^(-1) = E.
- n >= 0: by (LAM2) (with N = n+2), (E - (n+2)ψ) σ_0^(n+1) is the constant (n+2) Sum_j β_j^(n+1)/j, which vanishes by Z-infinity since n+1 >= 1, plus terms in Λ_{n+1} (the top leaf λ_{n+2} cancels because (E - (n+2)ψ) λ_{n+2} = 1). The Ŷ terms lie in Λ_{n-m+1}, and the m = n+1 term vanishes because σ_0^(-1) = 0. The E/t_c term of Ŷ is built into the prefactor E - (n+2)ψ.

So σ_1^(k) lies in Λ_{k+1} for k >= 0.

**Step 2 (σ_2).** Apply Lemma 4.3 with F = σ_1, using Step 1.
- For n <= -3 the coefficient [t_c^n] D_c σ_1 is 0: this uses σ_1^(-2) = 0, and it is what gives σ_2^(k) = 0 for k <= -2.
- For n = -2 it is E^2; for n = -1 it is a constant.
- For n >= 0 it lies in Λ^+_{n+1}, i.e. in Λ^+_{max(n+1, 0)}. A ψ-linear term could only come from (E - (n+2)ψ)·const(σ_1^(n+1)), which is zero: by Step 1 when n+1 >= 0, and because n + 2 = 0 for n = -2.

By (I2), using κ^(0) = μ^(0) = 0 and κ^(1) = -1: σ_2^(k) lies in Λ^+_k for k >= 0, and σ_2^(-1) = -E^2.

**Step 3 (constants of σ_2).** By Steps 1-2 each component of σ is bounded as ψ -> ∞, with limit equal to its constant part. In (K=1) the two arguments are ψ and -ψ; take the single limit ψ -> ∞. Using η^11 = E η^12, η^22 = 0 and σ_0 -> 0:

    -3 ε1 ε2 [E (E/t_c)^2 + 2 (E/t_c) σ_2^∞] = 3 ε1 ε2 E^3/t_c^2.

Hence σ_2^∞ = -E^2/t_c, and every constant of σ_2^(k), k >= 0, is 0. This limit is taken coefficientwise in t_c, and each coefficient is a rational function of ψ (Corollary 2.8).

**Step 4 ((c)).** With r_1 = 1/9 and v = -t_c/E - t_c^2/(18E) + O(t_c^3): σ_0^(1) = -(E/18) λ_1 + (E/9) λ_2. Then σ_1^(0) = (E - 2ψ) σ_0^(1) + (E^2/18) λ_1 = (E^2/9) λ_1. σ_2^(0) lies in Λ_0 = 0.

**Step 5 ((d)).** Coefficientwise in v: lim_{z->∞} S^[1] = S_∞'/E and lim_{z->∞} S^[2]/z = v S_∞''/E^2, S_∞ = Sum r_l (Ev)^{l+1}/(l+1)!. Lemma 4.1 then gives lim_{z->∞} P^_2/z = Δ_c/(C1i^2 v)·θ(C1i θ S_∞). By Z-infinity, S_∞ = -t_c, so C1i θ S_∞ = -C1cf/c1 = -1/f. Hence the formula, after multiplying by ρ and setting ψ = -z. ∎

Remark. A ψ-linear term Aψ in σ would add the ψ-constant -⟨A, P^^{(j)}(-w)⟩ to ρ Ehat_0j; this heuristic matches the off-shell controls, and is not load-bearing. Consistently deformed sections c1 = -wc(1 + Q Δ_c) were used as controls in round 1 and 3: the coefficient E Δ_c C1cf^2 θ(f)/t_c, indexed from k = 0 (value at k = -1 is 0), is matched exactly [E]. The earlier XXSHIFT controls shift XX without shifting c1, and are not a consistent section.

**Lemma G (gauge; R3:Lemma G).** [P] With the same c1 at both colours,

    η^11 P_1 P_1' + η^12 (P_1 P_2' + P_2 P_1')

depends only on (P_0, XX), because η^11 = E η^12. So every Ehat_ij is the same in the "hybrid gauge" (c1 := 1 in the prefactors but XX = θ log C1cf, as in edges_sing.pl; this is not the c1 = 1 polarisation) and in the true gauge. The map between the two gauges is z-free, analytic and invertible (C1cf(0) = -1), so H2 holds in either gauge. In particular the sign of c1 does not affect the edges. By s.1.4, the sign of c1 does not affect Dil^_i or Delta_i either.


## 5. The residue lemma and E_00

**Lemma R (residue lemma; R3:Lemma R; any rank, any constant symmetric pairing).** [P] Let f_0 = 1 and f_j = λ_j. Let P(ψ, ψ') = Sum_{i,l=0}^{K} C_il f_i(ψ) f_l(ψ') with C_il = C_li, and suppose κ + P(ψ, -ψ) = 0 for a constant κ. Then:
- κ = -C_00;
- E C_i0 + Sum_{l>=1} i C_il/(i+l) = 0 for every i >= 1;
- [κ + P(ψ, ψ')]/(ψ + ψ') = Sum_{i,l>=1} il/(E(i+l)) C_il λ_i(ψ) λ_l(ψ').

*Proof.* Use λ_i(ψ) - λ_i(-ψ') = i (ψ+ψ') λ_i(ψ) λ_{-i}(ψ') and λ_{-i} λ_l = (i λ_{-i} + l λ_l)/(E(i+l)). On the diagonal, P(ψ, -ψ) decomposes in the basis {1, λ_j : j != 0}, linearly independent because the poles are distinct. ∎

**Theorem E00 (R3:Theorem E00).** [P, using s.3.1 and H2 Steps 1-2] By (1.1) and K = 1, Ehat_00 = [Delta_0 - ⟨τ(ψ), τ(ψ')⟩]/(ψ + ψ'), τ(ψ) := P^(-ψ). Put 𝒴 := t_c σ and

    Ebar := t_c^2 ρ^2 Ehat_00 = [3 ε1 ε2 E^3 - ⟨𝒴, 𝒴'⟩]/(ψ + ψ').

Then:
- (i) Ebar = Sum_{n>=2} t_c^n Sum_{i,l>=1, i+l<=n} il/(E(i+l)) C^(n)_il λ_i λ_l', with C^(n)_il = -Sum_{k1+k2=n-2} Sum_{a,b} η^{ab} [λ_i]σ_a^(k1) [λ_l]σ_b^(k2).
- (ii) Ehat_00 lies in Q(ε1, ε2)[[ψ, ψ']][[t_c]]. In particular it is regular at the conifold.
- (iii) For all integers e and k >= 0: [t_c^k](ρ^e Ehat_00) lies in span{λ_i λ_l' : i, l >= 1, (i-1)+(l-1) <= k}, and it is symmetric.
- (iv) Ehat_00|_{t_c=0} = -(ε1 ε2 E^4/27) λ_1 λ_1'.
- (v) [t_c^k][ψ'^b](ρ^e Ehat_00) lies in Λ_{k+1}(ψ).

*Proof.* Apply Lemma R order by order in t_c, with κ = 3 ε1 ε2 E^3 δ_{n0}. Its hypothesis is (K=1). 𝒴^(0) = σ^(-1) = (0, E, -E^2) is constant. The constants of 𝒴^(a) feed only C_i0, C_0l, C_00, and these do not appear on the right of Lemma R. By H2, the leaf λ_i occurs in 𝒴^(a) only for i <= a, and λ_l occurs in 𝒴^(b) only for l <= b, with a + b = n. This gives (i), with k1 = a - 1, k2 = b - 1 >= 0. At n = 0, 1 there are no admissible pairs, giving (ii). (iii) and (v) follow by multiplying with the z-free unit ρ^{e-2}: at order n = k + 2, (i-1)+(l-1) <= n - 2 = k. For (iv): C_11 = -⟨b, b⟩ with b = (E, E^2/9, 0), ⟨b, b⟩ = 2 ε1 ε2 E^5/27, so the coefficient of t_c^2 in Ebar is (1/(2E)) (-2 ε1 ε2 E^5/27) λ_1 λ_1' = -(ε1 ε2 E^4/27) λ_1 λ_1', and ρ(0) = 1. ∎

Remarks.
- The hypotheses (H1)-(H3) of round 1 map as follows: H1 = Definition 1.3 (the identification with BS's rows, Theorem 2.7, is no longer needed here); H2 = Proposition H2; H3 = K = 1 (Proposition 2.10).
- Theorem E00 makes the rank-3 flow of the conifold s.7 unnecessary.
- Since the expansion (i) is explicit in t-independent leaves, no factorisation of the E_00 block is used; the issue of LSP Remark 13.1 does not arise.
- Only the Delta-normalised edge is treated here; vertex factors enter in s.7.


## 6. Regular colours and mixed edges (Theorem M)

Throughout, j, j' in {1, 2}, and v_{j,k} := [t_c^k] P^^{(j)}, a vector of formal power series in w (at z = -w). It is not a vector of rational functions: since v_j(0) != 0, the l-sum in S(v_j, -w) does not terminate (Corollary 2.8). Since Delta_j^{1/2} is analytic and non-zero at t_c = 0, Corollary 2.9(R) gives v_{j,k} = 0 for k < 0. (M-a)-(M-c) are taken [w^b]-wise, so they are not affected. Moreover, the component 2 of v_{j,0} is 0 identically: in Definition 1.3, P_2 has the factor Δ_c and the bracket is analytic at c1 = C1cf.

**Theorem M (R3:Theorem M).** [P, using Corollary 2.9(R), Proposition 2.10, H2]
- (M-a) ρ Ehat_0j and Ehat_0j lie in C[[t_c]][[ψ, w]].
- (M-b) Write σ^(k1) = c^(k1) + Sum_{l=1}^{k1+1} S_l^(k1) λ_l (c^(-1) = σ^(-1), S^(-1) = 0, and c^(k1) = 0 for k1 >= 0 by H2). Then

      [t_c^k](ρ Ehat_0j)(ψ, w) = -Sum_{k1=0}^{k} Sum_{l=1}^{k1+1} l ⟨S_l^(k1), v_{j,k-k1}(-w)⟩ λ_l(ψ)/(E + l w).

  At t_c = 0 this is rank one, λ_1(ψ) χ_j(w), with χ_j(w) = -⟨(E, E^2/9, 0), v_{j,0}(-w)⟩/(E + w) (from σ^(0), H2(c)).
- (M-c) For every integer e: [t_c^k][w^b](ρ^e Ehat_0j) lies in Λ_{k+1}(ψ), with no constant.
- (M-d) Ehat_j0(ψ, w) = Ehat_0j(w, ψ).
- (M-e) Ehat_jj', Dil^_j and Delta_j are analytic at the conifold; Delta_j(0) = ε1 ε2 E^3 ω^{2j}/(3(1-ω^j)^2) != 0; and Ehat_jj'(ψ, w) = Ehat_j'j(w, ψ).

*Proof.*

(M-a), (M-b). Put 𝒩 := ρ P_0j = -⟨σ(ψ), P^^{(j)}(-w)⟩. By (1.1), (ψ + w) ρ Ehat_0j = 𝒩(ψ, w) - 𝒩(-w, w). Every ψ-constant cancels in this difference: all of σ^(-1), and all constants at k >= 0 (which vanish anyway, H2). By (LAM3), the λ_l part contributes -(ψ + w) l ⟨S_l, v_{j,k2}(-w)⟩ λ_l(ψ) λ_l(-w), and λ_l(-w) = 1/(E + l w). Dividing by ψ + w gives (M-b). The coefficients for k <= -1 vanish because σ^(-1) is constant. By Corollary 2.8, each Ehat_0j(m, n) is a meromorphic germ with non-negative t_c-exponents, hence analytic.

(M-c). Take [w^b] in (M-b) and use the z-free unit ρ^{e-1}.

(M-d). This is the edge symmetry of s.1.4, from Proposition 2.10 at (i, j) = (0, j). That is ⟨P^^{(0)}(u), P^^{(j)}(-u)⟩ = 0, as an identity of Laurent expansions in t_c: by the drop-out, the pairing depends only on P_0 and each z-coefficient is rational in t. So P_0j(-w, w) = 0 = P_j0(-w, w), and η is symmetric.

(M-e). Each Ehat_jj'(m, n) is a polynomial in finitely many row coefficients, analytic by (R). Dil^_j uses component 0 only. The value of Delta_j uses 1 + 3t_j -> 1 - ω^j != 0, by s.1.1. The symmetry is the edge symmetry of s.1.4. ∎

Remarks.
- The regular-regular edge between two different colours has no conifold analogue; (M-e) covers it.
- Round 1 also found an exact s-recursion at the regular colours, S_k in s Q[s] of degree <= k+1, s = 1/(L + α), and a Lagrange identity for orthogonality. Both are superseded by the closed form (R2) and by Proposition 2.10, and are not used.
- The rank-one form at t_c = 0 was also verified directly in all entries [E, round-3 referee script R4], and the formula for χ_j reproduces the stored χ_1, χ_2 at (7,3) [E, round-4 critic 2]. Round 3's "conjectural rank-2 relation" and the inconsistent third entry of its stored vector are explained by a one-liner: the third entry multiplies a zero row. Indeed, v_{j,0} has component 2 equal to 0 by Definition 1.3, and the order-t_c^{-1} part of orthogonality with σ^(-1) = (0, E, -E^2) gives the same: ⟨σ^(-1), v_{j,0}⟩ = -3 ε1 ε2 E (v_{j,0})_2.

**Corollary 6.1 (edge inputs).** [P] Theorem LC_P2 holds, with no BS input (Proposition 2.10 is proved from Theorem A).

*Proof.* Table of sources:

    F1    Lemma 1.2(b), (M-e)
    F2    Corollary 2.9(F2), H2(a), (c)
    F3    Dil^_0: Lemma D0(iv).  Dil^_j: (M-e).
          Ehat_00: E00(ii).  Ehat_0j: (M-a).  Ehat_j0: (M-d).  Ehat_jj': (M-e).
    L1    Lemma D0(ii), (iii)
    L2    (M-c), (M-d); E00(v)
    L3    E00(iii), (iv)                                                         ∎

Values at t_c = 0, for orientation:

    factor     value at t_c = 0
    Dil^_0(d)  -E^{1-d}, i.e. q_t|_0 = q_0 = -Eψ/(E-ψ)
    Ehat_00    -(ε1 ε2 E^4/27) λ_1 λ_1'
    Ehat_0j    λ_1(ψ) χ_j(w), rank one, χ_j(w) = -⟨(E, E^2/9, 0), v_{j,0}(-w)⟩/(E + w)
    Ehat_jj'   analytic, no closed form given here
    Delta_0    3 ε1 ε2 E^3/(t_c ρ)^2 (double pole)
    Delta_j    ε1 ε2 E^3 ω^{2j}/(3(1-ω^j)^2) != 0


## 7. Hodge side and Theorem R_P2

Abbreviations: LSP = hodge_vertex_conifold/leaf_proof/leaf_structure_proof.md; PSM = hodge_vertex_conifold/proof_status_mixed_scale.md; CRN = hodge_vertex_conifold/conifold_regularity_note.md.

### 7.0 Dictionary conifold -> local P2

    conifold (LSP)                         local P2 (this text)
    rank 2                                 rank 3
    idempotent 1 (singular)                colour 0 (singular on the real branch; Lemma 2.1)
    idempotent 2 (regular)                 colours 1, 2
    t                                      t_c
    c = 1 + x                              E = ε1 + ε2
    lam_j = 1/(c - jψ)                     λ_j = 1/(E - jψ)
    q_0 = c - c^2 lam_1                    q_0 = E - E^2 λ_1 = -Eψ/(E-ψ)
    Delta_1 = -x c^3/(t ρ)^2               Delta_0 = 3 ε1 ε2 E^3/(t_c ρ)^2
    ρ = τ/(t alg)                          ρ = -E v_0/t_c
    r_{0,1}                                S(v_0, z) = component 0 of the normalised colour-0 row
    v_1 = (r_0, r_1)                       (P^_0, P^_1, P^_2)
    E_11; E_12, E_21; E_22                 Ehat_00; Ehat_0j, Ehat_j0; the four Ehat_jj'
    w = (c, -x, -1)                        w = (E, -ε1, -ε2) = ε2 (c, -x, -1) at x = ε1/ε2, using c = 1 + x (LSP App. A)

The CY relation c = 1 + x is what makes (E, -ε1, -ε2)/E = (1, τ, -1-τ) coincide with (c, -x, -1)/c, and what makes c = E at ε2 = 1. It is not a free substitution.

### 7.1 String and dilaton identities

For h >= 2 and a multiset R = {j_1, ..., j_r}, write ⟨f | R⟩ := Sum_m (1/m!) Int Omega_w f(ψ_0) prod_i λ_{j_i}(ψ_i) prod_k T(ψ_k).

**(H0) Dilaton homogeneity on 𝔇.** [P; this makes precise LSP s.13.1(i), which does not state a domain.] Let 𝔇 and 𝔡 be as in s.1.5. For q in 𝔇 and λ in 1 + t_c 𝐑 independent of ψ, one has λq in 𝔇 and

    d^n Fpot_h(λ q)[f_1, ..., f_n] = λ^{2-2h-n} d^n Fpot_h(q)[f_1, ..., f_n],

equivalently

    d^n Fpot_h(q)[f_1, ..., f_n] = λ^{2h-2} d^n Fpot_h(λ q)[λ f_1, ..., λ f_n].        (7.1)

Only λ = ρ is used. Omega_w is a rank-one CohFT at every colour, so nothing here depends on the rank.

*Proof.* Step 1 (base point). The dilaton equation for Omega_w, which holds because Lambda^vee is pulled back along forgetful maps, gives Int Omega_w ψ_0 prod g_i prod T = (2h - 2 + s + m) Int Omega_w prod g_i prod T on the spaces of s.1.5. Since q_0 = T - ψ, reindexing the T-sum gives

    𝔡^{s+1}(g_1, ..., g_s, q_0) = Sum_m [m - (2h - 2 + s + m)] A_m = (2 - 2h - s) 𝔡^s(g_1, ..., g_s),

where A_m := (1/m!) Int Omega_w prod g_i prod_{k<=m} T is the summand of 𝔡^s with m insertions of T. The m-dependence cancels. (At s = 0 this is E^2 Xi_h(0) - E^2 Xi_h(1) = (2 - 2h) Z[h,0], which is (L0_P2) + (C_P2).)

Step 2 (Euler identity on 𝔇). For q = q_0 + δ,

    d^{n+1} Fpot_h(q)[f, q] = Sum_r (1/r!) [𝔡^{n+r+1}(f, δ^r, q_0) + 𝔡^{n+r+1}(f, δ^{r+1})] = Sum_r (1/r!) (2 - 2h - n - r) 𝔡^{n+r}(f, δ^r) + Sum_r (1/r!) 𝔡^{n+r+1}(f, δ^{r+1}) = (2 - 2h - n) d^n Fpot_h(q)[f],

because the terms with the factor r telescope against the last sum.

Step 3 (integration). Put μ := λ - 1, in t_c 𝐑. Let ϑ be a formal variable, q_ϑ := (1 + ϑμ) q, and g(ϑ) := d^n Fpot_h(q_ϑ)[f]. Each t_c-coefficient of g is polynomial in ϑ. The Taylor definition gives dg/dϑ = μ d^{n+1} Fpot_h(q_ϑ)[f, q], and Step 2 at q_ϑ gives (1 + ϑμ) dg/dϑ = μ (2 - 2h - n) g. The solution with value g(0) at ϑ = 0 is unique, so g(ϑ) = (1 + ϑμ)^{2-2h-n} g(0). Put ϑ = 1. ∎

**(S_P2) String identity.** [P] Xi_h(1, R) = (J_R/E^2) Xi_h(R). Iterating, Xi_h(1^p, R) = prod_{i=0}^{p-1} (J_R + i) · E^{-2p} Xi_h(R).

*Proof.* By (LAM4), Xi(1, R) = E^{-1} [⟨1 | R⟩ + ⟨ψ λ_1 | R⟩]. The string equation (valid since h >= 2 and Omega_w is pulled back) gives S λ_j = (λ_j - λ_j(0))/ψ = (j/E) λ_j on probes, and S T = T/ψ = -ψ λ_1 on dilaton leaves. Summed over m, ⟨1 | R⟩ = (J_R/E) Xi(R) - ⟨ψ λ_1 | R⟩, and the second term cancels. ∎

**(C_P2) Theorem C.** [P] Xi_h(1^p) = 0 for h >= 2, p >= 1, i.e. Z[h, n] = 0 for n >= 1. *Proof.* (S_P2) with R empty: J = 0 kills every p >= 1. It uses only the string equation and the leaf shapes, no rank. (This is CRN s.3a with c -> E.) ∎

**(L0_P2) The λ_0 identity.** [P] Xi_h(0, R) = ((D_R - 2h + 2)/E^2) Xi_h(R), with λ_0 = 1/E.

*Proof.* ψ λ_1 = (ψ - T)/E (from T = E - E^2 λ_1 + ψ). The dilaton equation gives ⟨ψ | R⟩ - ⟨T | R⟩ = (2h - 2 + r) Xi(R): after reindexing the T-sum, the 1/m! cancels the m-dependence. So ⟨1 | R⟩ = ((J_R - r - 2h + 2)/E) Xi(R). Multiply by λ_0 = 1/E. ∎

Evidence: checks.py and ref2.py [E, Appendix]. These run on the trivial CohFT at genus 0 and 1 and are controls of the algebra, not of the Hodge statement; the proofs above are complete.

### 7.2 Why λ_0 would be fatal and does not occur

A constant component at an insertion of order t_c^k cannot be stripped by (S_P2). Its D-value is -1: it lowers D instead of being absorbed by the counting. At R' empty it leaves ((2 - 2h)/E^2) Z[h, 0], a gap multiple, at order 2 - 2h + k. So a λ_0 component would be a genuine obstruction. It does not occur:
- H2(b): σ^(k), k >= 0, has no constants and no ψ-polynomial part;
- E00(i): only λ_i λ_l with i, l >= 1 appear;
- (M-c): no constant;
- Lemma D0(iii): no constant for k >= 1.

The constant E in q_0 is the base point of the Taylor expansion, not an insertion. σ^(-1) is a constant, but it only ever feeds C_i0, C_0l, C_00 in Lemma R and it cancels in (M-a); it never reaches a vertex.

### 7.3 Range of Lemma A

**Lemma 7.1.** [P] Let Γ be a connected stable graph of genus g, without legs, with at least one edge. Then every g_v <= g - 1. In particular Theorem R_P2 uses M_irr(h) only for 2 <= h <= g - 1, and at g = 2 there is no rescaled vertex at all.

*Proof.* If h^1(Γ) >= 1 this is immediate from g = h^1(Γ) + Sum g_v. If Γ is a tree with >= 2 vertices, it has at least 2 leaf vertices of valence 1. Each has genus >= 1 by stability. Hence no vertex carries all of g. ∎

The stability count concerns half-edges only: for G_g there are no external markings, and the "legs" of step (ii) below appear only after the graph has been fixed. The decorations agree with the code: hv_lib3.m vertexOptionsH (lines 119-130) has dilaton degrees >= 2 and dim >= 2 nd, so a genus-0 vertex with m = 0 needs at least 3 half-edges. (Round-1 note N6: the single-vertex graph does need h = g.)

### 7.4 Transfer of Lemma A (Lemma A_P2)

**Lemma A_P2.**

(a) *Identification.* [P] Xi_h(R) of s.1.5 is the conifold object of Lemma A: same normalised class, same leaves, same dilaton point, with c -> E. In the normalisation c = 1 of Lemma A (PSM conventions; lemmaA/round2/brief.md s.1: w = (1, τ, -1-τ), Omega = Lambda^vee(1) Lambda^vee(τ) Lambda^vee(-1-τ)/(-σ_τ)^h, lam_j = 1/(1 - jψ), T = -ψ^2/(1-ψ), σ_τ := τ(1+τ)), one has (E, -ε1, -ε2) = E (1, τ, -1-τ), τ := -ε1/E, -1-τ = -ε2/E, σ_τ = -ε1 ε2/E^2. And

    Xi_h(R)(ε1, ε2) = E^{-(3h-3+2r)} Xi_h^{(c=1)}(R)|_{τ = -ε1/E}.       (7.2)

*Proof of (7.2).* Weight count. Write Omega_w = Sum_i Omega_i with Omega_i of cohomological degree i, homogeneous of degree -i in w (this uses the normalisation by w_a^h; with the unnormalised product the exponent would be 3 - 2r). The degree-a part of a leaf 1/(E - jψ) scales as E^{-1-a}, and of a dilaton leaf -ψ^2/(E - ψ) as E^{1-a}. On Mbar_{h,r+m} the total degree is 3h - 3 + r + m. Hence the factor is E^{-(3h-3+r+m) - r + m} = E^{-(3h-3+2r)}, independent of m. ∎

(b) *Vanishing.* [C: hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md, Theorem L and Corollaries L1-L3, s.7.5 there, for all h >= 2; there P modulo (MV) and textbook facts; status below] For every h >= 2 and every non-empty multiset R of integers >= 2 with D(R) <= 2h - 3: Xi_h(R) = 0 in Q(ε1, ε2), and at every specialisation with ε1 ε2 E != 0, the diagonal ε1 = ε2 (τ = -1/2) included. Combined with (S_P2) and (C_P2), the same vanishing holds for every non-empty R with parts >= 1 in the range. R must not contain a part 0: Xi(0) = ((2-2h)/E^2) Z[h,0], which is a non-zero element of Q(ε1, ε2) by Theorem VI.

*Proof.* Theorem L of LEMMA_A_PROOF.md says Xi_h^{(c=1)}(R) = 0 in Q(τ) for such R, with R non-empty. (R3's own statement omits "non-empty"; the correction is round2_status.md s.1, "precise form".) By its Corollary L1 (= PSM s.1.1(e)), Xi_h^{(c=1)}(R; τ) = Sum_{k=0}^{h} (-1/σ_τ)^k Xi^(k)_h(R), with τ-independent Xi^(k)_h. So the identity specialises at any τ with σ_τ != 0, and (7.2), which is its Corollary L3, gives (b). ∎

**Status of the proof of Lemma A.** The written proof is hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md, consolidating the round-5 reports ThmL and MVmatch (journal wf_9d7fd94e-8d7).
- **Proof of record:** R3 proof (ii), all u_k-linear members P_k, rewritten uniformly (parts >= 1, no empty-profile table), with the bilinear equations P_k derived for all k in s.4 of that file. It cites no KP or 2D Toda theorem, and uses no cut-and-join, no Route A Theorem B or Cor. E, no TR, and no computer input. This answers round-4 critic 1 S2.
- **Matching with Xi_h:** report MVmatch (LEMMA_A_PROOF.md Proposition 2.1 and s.3; all-profile form Proposition 8.1) derives the matching of the MV series with Xi_h for all genera and all profiles, every constant derived, with no data comparison. This answers round-4 critic 1 S1, i.e. item (iii) of the round-4 list. It does not by itself discharge the other inputs, which is why the proof of record above is listed separately.
- **Inputs that remain literature:** (MV), the Marino-Vafa formula, LLZ eq. (1) (eqs. (1)-(3) of arXiv:math/0306434v3; published JDG numbering not confirmed; "LLZ Thm 1" in earlier notes is a misnumbering, since LLZ's Theorems 1, 2 are the cut-and-join equations), with OP Theorem 1 as a second published proof (with the caveat on its printed 1/|Aut μ|, LEMMA_A_PROOF.md s.2.2); and textbook facts (string and dilaton equations for pulled-back classes, the genus-0 multinomial formula, Frobenius formula, Murnaghan-Nakayama, the hook-length form of the principal specialisation, Lagrange inversion). What was and was not confirmed is LEMMA_A_PROOF.md s.9.
- **Refereeing:** each of the two source reports had two referee passes (logic/literature and independent computation), all "sound with minor fixes", no major item. The assembled LEMMA_A_PROOF.md has not had its own referee (its open point O1). That file's status table (s.10) takes precedence over this summary. Finite-order runs (its s.12) are evidence only.
- **Not provided by that file:** the identification of the Taylor coefficients of the graph sum with Xi_h(R), which is s.7.5(ii) here, and the empty profile, which is Theorem VI here. Its open point O5 and its s.8.5, written in parallel with this edit, still list (VI) and lem:modrel as resting on R2 Cor. 7 and on critic 1 S4; s.7.8 and Proposition 2.12 of this text supersede that.

Parameter dependence, audited by reading:
- The induction lives in Q[τ][[x]] and Q(τ). Its divisions are by polynomials not identically zero: Π_2 = 2τ + 1 (which does vanish on the diagonal τ = -1/2) and Π_k. This is harmless, because **every inductive hypothesis and conclusion is an identity in Q(τ), never used at a specialised τ**. The specialisation happens only at the very end, through the Laurent polynomial in σ_τ above.
- It requires no sign, rationality or integrality condition. Only σ_τ != 0 and E != 0 are needed, which on T_F are the standing assumptions.
- Backup: round 2 route R2 (Lemma D proved at integer τ and extended by polynomiality; TR on the framed-C^3 curve C_τ, Hodge side, allowed; C_τ is not the local-P2 mirror curve). Its ABDKS/Eynard normalisation is not derived for all g (round2_status.md s.5.1 item 6), so it is a cross-check, not a proof of record.
- **R1 is not admissible** under the all-orders ruling Q5: its base h <= 7 and its zero count are computer-assisted. The computations at h = 3, 4 (PSM s.1.9) are likewise unnecessary.

(c) *Range.* [P] Theorem R_P2 uses (b) for 2 <= h <= g - 1 (Lemma 7.1). Corollary SV_P2 also uses (b) at h = g. Both are covered.

The fact that makes this transfer possible is not automatic and is checked in s.1.6: the leaf scale E produced by F2 equals the first Hodge weight w_1 = E.

### 7.5 Proof of Theorem R_P2

Fix Γ and a colouring c. Call a colour-0 vertex with g_v >= 2 **rescaled**. For a graph with an edge every vertex has at least one half-edge. Fix ρ as in s.1.2. Dilaton decorations at every vertex are summed first, as in (1.2). Everything below is multilinear in the half-edge slots: for an edge that does not factorise (a loop, or an edge between two rescaled vertices), the statements are applied term by term in the expansion of E00(i) or E00(iii), and they are read coefficientwise in the edge variables of the far end.

**(i) Homogeneity.** [P] At a rescaled vertex v of genus h with half-edge series f_1, ..., f_n, by (1.2), (7.1) with λ = ρ (q_t lies in 𝔇 by Lemma D0(iv), and ρ q_t lies in 𝔇 although its ψ-coefficient is -ρ), and Lemma 1.2(b), the vertex factor is

    Delta_0^{h-1} d^n Fpot_h(q_t)[f_1, ..., f_n] = (3 ε1 ε2 E^3/t_c^2)^{h-1} d^n Fpot_h(ρ q_t)[ρ f_1, ..., ρ f_n].

Only two properties of ρ are used: it is an analytic z-free unit, and ρ^2 Delta_0 is a constant times t_c^{-2}.

**(ii) Leaf expansion.** [P] By Lemma D0, ρ q_t = q_0 + Sum_{k>=1} t_c^k q^(k) with q^(k) in Λ_{k+1}. Each half-edge at v receives one factor ρ:
- *edge to a regular colour:* [t_c^k][w^b](ρ Ehat_0j) lies in Λ_{k+1} with no constant, for both orientations: (M-c), (M-d);
- *edge to an unrescaled colour-0 vertex (genus 0 or 1):* ρ Ehat_00; [t_c^k][ψ'^b] of it lies in Λ_{k+1}: E00(v) with e = 1;
- *edge between two rescaled vertices, or self-loop:* ρ^2 Ehat_00, whose order-t_c^k coefficient lies in span{λ_i λ_l' : (i-1)+(l-1) <= k}: E00(iii). The leaf degree bound is split between the two ends, or assigned wholly to v for a loop.

At an unrescaled vertex only regularity is needed:
- Delta_0^{-1} = (t_c ρ)^2/(3 ε1 ε2 E^3) = O(t_c^2) at colour-0 genus 0;
- Delta_0^0 = 1 at colour-0 genus 1;
- Delta_j^{h-1} is analytic for every h, including h = 0, because Delta_j(0) != 0 (M-e);
- their dilaton and edge factors are analytic (Lemma D0(iv), (M-e), E00(ii));
- the unit ρ^{-1} does not raise leaf degree.

By the definition of Fpot_h on 𝔇 (s.1.5), Fpot_h(q_0 + δ) with δ = Sum_{k>=1} t_c^k q^(k) is its Taylor series at q_0. The expansion is finite at each order, because δ = O(t_c) and each Xi is a finite sum. So the coefficient of t_c^{2-2h+k} of the rescaled vertex factor is a finite combination of Xi_h(R) with all parts j >= 1 and D(R) <= k. There are no λ_0 terms and no ψ-polynomial remainders (s.7.2). Each half-edge contributes at least one leaf, so r >= n_v >= 1.

All leaves at order t_c^k are λ_j with 1 <= j <= k + 1: only lattice points of M enter. Since H2 and E00 are explicit expansions in t-independent leaves, no factorisation is used (LSP Remark 13.1 does not arise).

**(iii) Vanishing.** [P mod Lemma A_P2(b)] Strip the parts equal to 1 by (S_P2): Xi_h(1^p, R') = prod_{i<p} (J_{R'} + i) E^{-2p} Xi_h(R'), where R' has parts >= 2.
- If R' is empty, then p = r >= 1 and J_{R'} = 0 kills the term (C_P2).
- If R' is non-empty, D(R') = D(R) <= k. If k <= 2h - 3 the term vanishes by M_irr(h), with h <= g - 1 by Lemma 7.1.

So a non-zero term needs D_v >= 2 h_v - 2 at every rescaled vertex, where D_v is the total leaf degree at v.

**(iv) Counting.** [P] Consider a non-zero term in the expansion of T(Γ, c). Each factor contributes a t_c-order, and its leaf degree at any rescaled end is bounded by that order (with an edge between two rescaled vertices split, and with a loop assigned wholly to v). So the t_c-order of the term is at least

    Sum_{v rescaled} (2 - 2 h_v + D_v) + 2 #{colour-0 vertices of genus 0} + 0,

where the final 0 accounts for the regular-colour vertices, the genus-1 colour-0 vertices, and every edge/dilaton order not allocated to a D_v, all of which contribute >= 0. By (iii) each summand is >= 0. Hence the order is >= 0. ∎

**Status of the conclusion.** Theorem R_P2 holds for the decoration sum as defined (s.1.4), per colouring, for every genus, modulo Lemma A_P2(b) for 2 <= h <= g-1. Genus 2 uses no Lemma A at all. No BS input enters (s.0.3). Identification with GW needs (C1)-(C5) and Proposition 13.2_P2, which remains open.

### 7.6 The single-vertex graph; Corollaries SV_P2 and G

**Theorem VI (vertex identity (VI)).** [P mod LLZ eqs. (1), (5), the boson-fermion correspondence, string/dilaton; proof in s.7.8] For g >= 2,

    (-ε1 ε2 E^3)^{g-1} Z[g, 0] = β_{2g-2}(ε1, ε2),     β_n := (n-1)! [t^n] 1/(4 sinh(ε1 t/2) sinh(ε2 t/2)).

At ε2 = 1, ε1 = x, with c = 1 + x = E, this reads (-x c^3)^{g-1} Z[g, 0] = gap_g, where gap_g = (2g-3)! [t^{2g}] t^2/(4 sinh(xt/2) sinh(t/2)) = (2g-3)! P_g(x,1)/x. This is conifold F4 at n = 0.

Z[g, 0] is a pure Hodge integral: it depends only on Omega_w and the scale E. The proof in s.7.8 uses no TR, no Eynard/ABDKS, no Lemma A, no Teleman and no finite-order input.

Cross-checks, not proof steps:
- (B-a) Lemma A round 2, route R2 Corollary 7 (gap formula, lemmaA/round2/round2_status.md s.1.2), by Hodge-side TR. It is equivalent to Theorem VI by the short conversion of s.7.8, Part A. Its ABDKS/Eynard normalisation was not derived for all g, so it is demoted to a cross-check.
- (B-c) conifold Proposition 13.2 + BS prop:ResConifold (main.tex line 1708; "Prop. 3.1" of arXiv:2410.00118) (H) + Corollary R with Lemma A(<= g-1). This is TR-free, but the Teleman input of conifold Prop 13.2 was read through summaries only (conifold critic gap list), and it keeps a BS citation. So it is a consistency remark. See also BS prop:cfgraphsum (statement at 5114, eq:betagiv 5122, label 6366), which is related and which has not been reconciled with the single-vertex value Z[g,0].

The recorded proof is s.7.8.

**Proof of Corollary SV_P2.** [P mod M_irr(h), h <= g, and Theorem VI]
- *Colours 1, 2.* Delta_j^{g-1} and Dil^_j are analytic: the term is regular.
- *Colour 0.* By (i)-(iii) with no half-edges, the vertex factor is (3 ε1 ε2 E^3/t_c^2)^{g-1} Fpot_g(ρ q_t). Its k-th Taylor term is a combination of Xi_g(R) with R empty only at k = 0, and R non-empty with parts >= 1, D(R) <= k otherwise. For 1 <= k <= 2g - 3 these vanish by Lemma A_P2(b) at h = g (with R' empty killed by C_P2). So the polar part is exactly (3 ε1 ε2 E^3)^{g-1} Fpot_g(q_0) t_c^{2-2g}.
- Since 3 ε1 ε2 E^3 = (-3)(-ε1 ε2 E^3), (VI) gives (-3)^{g-1} β_{2g-2} t_c^{2-2g}.
- Theorem R_P2 removes all other polar terms from G_g, since (M_irr(h), h <= g) covers its hypothesis. ∎

The (-3)^{g-1} ratio is thus explained: same vertex value, Delta-prefactor 3 ε1 ε2 E^3 against -x c^3.

Consistency with BS [C, not a step]:
- BS thm:cfasym (5301): leading polar coefficient (-3)^{g-1} β_{2g-2}/t_cf^{2g-2}.
- BS lem:cfpole (5407): pole order <= 2g - 2.
- BS rmk:cfgap (5422): the same formula with remainder O(1) is conjectural, verified at g = 2. At g = 2 our formula reads -(7ε1^4 + 10ε1^2ε2^2 + 7ε2^4)/(1920 ε1 ε2 t_c^2), exactly BS's display.

**Corollary G.** [P mod Corollary SV_P2 and Proposition 13.2_P2 [O]] If G_g = GW^cf_g := Cf(GW_g), then by Lemma 1.2(c), GW^cf_g = (-3)^{g-1} β_{2g-2}/t_cf^{2g-2} + O(1). This is rmk:cfgap. It is stated conditionally, and the condition is open (s.8.3).

### 7.7 Audit of the conifold proof for rank-2 or conifold-specific assumptions

- **(a) "Every non-singular vertex is analytic".** There are now two regular colours. Analyticity needs Delta_j(0) != 0 at both j, because regular vertices of genus 0 carry Delta_j^{-1}. It also needs Ehat_12, Ehat_21 between different regular colours (no conifold analogue). Both are in (M-e).
- **(b) Colouring sums.** None are used; every step is per colouring. Each T(Γ, c) contains no square root (s.1.4).
- **(c) Per-colouring statement.** Only the label of the singular colour depends on the branch. Near y = -1/27 the roots t_i are analytic and distinct (Lemma 2.1), so labels are locally constant. Continuation is along the real branch y in (-1/27, 0) (round-1 N8, journal wf_e928f771-9e3).
- **(d) Exclusion of the single-vertex graph.** It enters only through r >= 1 in (iii) and through h <= g-1 (Lemma 7.1).
- **(e) Rank.** Theorem R uses no 2x2 structure. Rank enters only inside the factor theorems: (1.1) and Lemma R are rank-free; Lemma G uses η^11 = E η^12; the orientation clause comes from (M-d).
- **(f) ε1, ε2 indeterminate or complex.** Theorem A and the proof of Proposition 2.10 use only E != 0 (and a square root of 3 ε1 ε2 E in S0), so they work verbatim over C(t) for every complex ε with ε1 ε2 E != 0. The statements of s.3-s.6 follow algebraically. Lemma A and Theorem VI hold identically in τ and specialise, through the Laurent polynomial of s.7.4(b) and through (7.2). But the non-vanishing of Z[h,0] does not specialise: at h = 2, Z_2(τ) = -(24 + 28s + 7s^2)/5760 with s = 1/σ_τ, by Theorem VI, and this vanishes at two real values σ_τ < -1/4, i.e. at complex τ. No proof step uses the non-vanishing.
- **(g) Polarisation.** Only XX enters the edges (Lemma G). Dil^ and Delta do not depend on c1 (s.1.4). On-shell is essential through H2; but the Hodge-side reduction and the counting do not depend on it.
- **(h) Lf factors.** They occur only with marked points, so not in G_g.
- **(i) Theorem R as stated in LSP/PSM** assumes M_irr(h) up to h = g; for edge graphs only h <= g - 1 is needed (Lemma 7.1). The documentation of the conifold should be corrected accordingly (round2_status.md A4).

### 7.8 Proof of Theorem VI

Normalisation c = 1 as in Lemma A_P2(a): τ an indeterminate, σ := σ_τ = τ(1+τ), σ' := dσ/dτ = 1 + 2τ, Omega = Lambda^vee(1) Lambda^vee(τ) Lambda^vee(-1-τ)/(-σ)^h, T = -ψ^2/(1-ψ), a_1 := -τ, a_2 := 1 + τ. So a_1 a_2 = -σ, a_1^2 + a_2^2 = 1 + 2σ, τ a_2 = σ. Put Z_h(τ) := Sum_m (1/m!) ⟨Omega prod T(ψ_k)⟩_{h,m}, which is Fpot_h(q_0) at c = 1.

**Theorem VI_τ.** For every h >= 2, identically in τ: σ^{h-1} Z_h(τ) = β_{2h-2}(a_1, a_2).

**Inputs.** LLZ eq. (1)-(3) and Theorem 1 eq. (5) (arXiv:math/0306434v3, read; no JDG numbering confirmed); the boson-fermion correspondence (Date-Jimbo-Kashiwara-Miwa; Miwa-Jimbo-Date, Solitons; no theorem number confirmed or asserted, and the derivation of the bilinear identity for this g is written out in Step 2(b), so the citation carries only the vertex-operator form of ψ, ψ^*); the principal specialisation s_ν(q^ρ) = q^{κ_ν/4}/prod_box [h] (Macdonald I.3 Ex. 2, recalled, example number not confirmed; LLZ eq. (13) states it; the exponent was re-derived by hand); string and dilaton equations for pulled-back classes; (7.2). The proof combines the MV formula, the first mKP equation (KP/Toda is allowed by the rulings) and the cut-and-join eigenvalue on the p_1-slice, extracts the lowest order at the dilaton point, and uses the τ-symmetry of Omega to fix the integration constants.

**Part A. Homogenisation.** [P] Theorem VI_τ implies Theorem VI. (7.2) at r = 0 gives Z[h, 0] = E^{-(3h-3)} Z_h at τ = -ε1/E, where a_i = ε_i/E. Further -ε1 ε2 E^3 = σ E^5 and β_n(λa, λb) = λ^n β_n(a, b). So (-ε1 ε2 E^3)^{h-1} Z[h,0] = σ^{h-1} E^{5h-5-3h+3} Z_h = E^{2h-2} β_{2h-2}(a_1, a_2) = β_{2h-2}(ε1, ε2). ∎

Conversion to R2 Cor 7 (for the record; not used). With S_R(y) = (y/2)/sinh(y/2), one has S_R(a l) S_R(b l) = a b l^2/(4 sinh(a l/2) sinh(b l/2)). So R2 Cor 7, σ^h Int K Omega = -(2h-3)! [l^{2h}] S_R(a_1 l) S_R(a_2 l), reads σ^{h-1} Int K Omega = β_{2h-2}(a_1, a_2). And PSM s.1.1(a) at r = 0 gives Sum_m (1/m!) π_* prod T(ψ_k) = exp(-Sum κ_k/k) = K, because 1 - T(z)/z = 1/(1-z). Hence Int K Omega = Z_h, and Theorem VI_τ is R2 Cor 7. This replaces the [E] gap_identity_check.pl claim of round 4.

**Part B. Proof of Theorem VI_τ.** Notation: M_h(1^a) := ⟨Omega prod_{i=1}^a 1/(1-ψ_i)⟩_{h,a}, G_h(y) := Sum_a M_h(1^a) y^a/a!. In genus 1 the sum starts at a = 1. M_0(1^a) := a^{a-3}, the LLZ unstable convention.

*Step 1 (the MV formula in real form, all g and all μ).* [C LLZ eqs. (1)-(3); the conversion is P] LLZ state: Sum_{μ,g} λ^{2g-2+l} C_{g,μ} p_μ = log Sum_ν e^{i(τ+1/2) κ_ν λ/2} V_ν(λ) s_ν(p), with C_{g,μ} = -(i^{|μ|+l}/|Aut μ|) σ^{l-1} Π_μ I_g(μ), I_g(μ) = Int Lambda^vee(1) Lambda^vee(-τ-1) Lambda^vee(τ)/prod (1 - μ_i ψ_i), where l = l(μ), Π_μ = prod_i Π_{μ_i}, Π_m = prod_{a=1}^{m-1}(mτ + a)/(m-1)!, and V_ν is the product of LLZ eq. (3),

    V_ν = prod_{a<b} (sin[(ν_a - ν_b + b - a)λ/2]/sin[(b - a)λ/2]) prod_i prod_{v=1}^{ν_i} 1/(2 sin[(v - i + l)λ/2]).

The hook form follows: with [n] := 2 sin(nλ/2) and l_i := ν_i + l - i, one has prod_box [h] = prod_i [l_i]!/prod_{a<b} [l_a - l_b] and prod_i [l - i]! = prod_{a<b} [b - a], so V_ν = prod_box 1/(2 sin(hλ/2)). Put x = iλ. Then 2 sin(hλ/2) = -i·2 sinh(hx/2), so V_ν = i^{|ν|} prod_box (2 sinh(hx/2))^{-1}, and i^{|ν|} s_ν(p) = s_ν(i^k p_k). On the left, λ^{2g-2+l}(-i^{|μ|+l}) = (-1)^g i^{|μ|} x^{2g-2+l}. Cancelling i^{|μ|} and using I_g = (-σ)^g M_g gives, with 𝖹(p) := Sum_ν e^{(τ+1/2) κ_ν x/2} prod_box (2 sinh(hx/2))^{-1} s_ν(p) and 𝖥 := log 𝖹:

    [x^{2g-2+l} p_μ] 𝖥 = σ^{g+l-1} Π_μ M_g(μ)/|Aut μ|.   (B1)

This derivation is uniform in g and μ. By the principal specialisation, 𝖹 = Sum_ν e^{τ x f_2(ν)} s_ν(q^ρ) s_ν(p), with f_2 = κ/2, q = e^x, [h] = q^{h/2} - q^{-h/2}. On the slice p = (p_1, 0, 0, ...), set p_1 = y/(σx). Then (B1) gives

    𝖥 = Sum_{h>=0} σ^{h-1} x^{2h-2} (G_h(y) - G_h(0))   in Q(τ)[[y]]((x)).   (B2)

There is no a = 0 term on the Schur side, since 𝖹(0) = 1. That constant plays no role below.

*Step 2 (a closed equation on the p_1 slice).*

(a) Cut-and-join. [C: LLZ Theorem 1 eq. (5), which LLZ attribute to Zhou; the rescaling p_k -> i^k p_k does not affect it, since the operator preserves degree] Let 𝒞 := (1/2) Sum_{i,j} (i j p_{i+j} ∂_i ∂_j + (i+j) p_i p_j ∂_{i+j}), ∂_i = ∂/∂p_i. Then 𝒞 s_ν = f_2(ν) s_ν, so ∂_τ 𝖹 = x 𝒞 𝖹. On the slice, the first sum of the logarithmic form vanishes (it carries p_{i+j} with i + j >= 2), and the second gives

    ∂_τ 𝖥|_p = x p_1^2 ∂_2 𝖥.

(b) First mKP. [P, with C for the boson-fermion dictionary] With t_k = p_k/k and g = e^{τ x F_2} Γ_-(p'), p'_k = 1/[k], put τ_n(t) := ⟨n| Γ_+(t) g |n⟩ (half-integer modes, |n⟩ filled below n). The charge-n eigenvalue of F_2 on |ν, n⟩ is f_2(ν) + n|ν| + C(n), so τ_0 = 𝖹 and τ_1 = e^{τ x/8} 𝖹(e^{kβ} p_k) with β := τ x. Since (A⊗1 + 1⊗A) commutes with Sum_k ψ_k ⊗ ψ_k^* for bilinear A, and Sum_k ψ_k|1⟩ ⊗ ψ_k^*|0⟩ = 0 (ψ_k|1⟩ = 0 for k < 1, ψ_k^*|0⟩ = 0 for k > 0), we get Sum_k ψ_k g|1⟩ ⊗ ψ_k^* g|0⟩ = 0. The vertex-operator form of ψ, ψ^* then gives Res_z z e^{ξ(t-t',z)} τ_1(t - [1/z]) τ_0(t' + [1/z]) dz = 0. Setting t = a + b, t' = a - b and taking the coefficient at b = 0 (the z^{-2} coefficient is (1/2)(D_1^2 - D_2) τ_1 · τ_0, D_k the Hirota derivatives) gives

    (D_1^2 - D_2) τ_1 · τ_0 = 0.

The form is derived, not fitted. In logs, with 𝖥_1(p) := 𝖥(e^{kβ} p_k):

    2 ∂_2(𝖥 - 𝖥_1) + ∂_1^2(𝖥 + 𝖥_1) + (∂_1(𝖥 - 𝖥_1))^2 = 0.

(c) Elimination of ∂_2. On the slice, ∂_2 𝖥_1 = e^{2β} (∂_2 𝖥)(e^β p_1), and (a) gives ∂_2 𝖥(P) = ∂_τ 𝖥(P)/(x P^2), so the e^{2β} cancels. With D := y d/dy = p_1 ∂_1, multiplying by p_1^2 gives

    (2/x)(1 - e^{βD}) 𝖥_τ + (D^2 - D)(1 + e^{βD}) 𝖥 + ((1 - e^{βD}) D𝖥)^2 = 0,   (*)

where 𝖥_τ := ∂_τ|_p 𝖥 = (∂_τ|_y + (σ'/σ) D) 𝖥, because y = σ x p_1, and where ∂_τ acts before e^{βD} and never differentiates β. (*) is an exact identity in Q(τ)[[y]]((x)). It contains 𝖥 only through D𝖥 and through D-derivatives of 𝖥_τ, so constants in y are irrelevant.

*Step 3 (the dilaton point).* [P, given string and dilaton] Lambda classes pull back under forgetful maps. For h >= 2 put FF(c, p, r) := Sum_m (1/m!) ⟨Omega prod (c + rψ + p/(1-ψ))(ψ_i)⟩_{h,m}. String and dilaton read (1-r) ∂_c FF = p ∂_p FF and (1-r) ∂_r FF = (2h - 2 + c ∂_c + p ∂_p) FF. These determine the c- and r-expansions recursively from FF(0, p, 0) = G_h(p), so FF(c, p, r) = (1-r)^{2-2h} G_h((p/(1-r)) e^{c/(1-r)}). Since ψ^2/(1-ψ) = 1/(1-ψ) - 1 - ψ, taking (c, p, r) = (-s, s, -s) gives Q_h(s) := Sum_k (s^k/k!) ⟨Omega (ψ^2/(1-ψ))^k⟩, and

    G_h(y) = w^{2h-2} Q_h(w - 1),   w := 1 + s = 1/(1-u),   y = u e^{-u}.

Q_h has degree <= 3h - 3 and Q_h(-1) = Z_h. So G_h = w^{2h-2}(Z_h + O(w)), a polynomial in w. Also D = (w-1) w^2 d/dw; with v := 1/w this is D = (1 - 1/v) d/dv.

Genus 0: D^2 G_0 = Sum a^{a-1} y^a/a! = u = 1 - v, hence G_0 = 5/12 - v^2/4 - v^3/6.

Genus 1: with a marked leg 1/(1-ψ_0), the string and dilaton equations become (1-r) ∂_c FF = (p ∂_p + 1) FF and (1-r) ∂_r FF = (1 + c ∂_c + p ∂_p) FF. So FF = (1-r)^{-1} e^{c/(1-r)} H((p/(1-r)) e^{c/(1-r)}), with H(y) := Sum_a M_1(1^{a+1}) y^a/a!. At (-s, s, -s) this gives H = w e^u Q_1(s; (1)), hence DG_1 = yH = u w Q_1(s; (1)), where Q_1(s; (1)) = Int_{1,1}(ψ + Omega_1) + s ⟨τ_2 τ_0⟩_1 = (s - 1/σ)/24 (in genus 1, Omega = 1 - ((1+σ)/σ) λ_1). So DG_1 = (w-1)(w-1-1/σ)/24, i.e. G_1 = (1/24)(-log v + (1 + 1/σ) v) + const.

Every x-coefficient of D𝖥 is a Laurent polynomial in w, and Laurent polynomials embed in Q(τ)[[y]] because w = 1 + y + ... is non-constant. So (*) holds coefficientwise in the monomials x^a w^b.

*Step 4 (grading).* [P] Map x^a w^b to ε^{b-a} N^{-b}, with ε := 1/x and v = xN. Each coefficient has a unique preimage, and N-powers are bounded above at fixed ε-order. This monomial substitution intertwines D with ε(1 - ε/N) d/dN, and βD with τ(1 - ε/N) d/dN, a derivation. So e^{βD} f(N) = f(Φ(τ)), where dΦ/dϑ = 1 - ε/Φ, Φ(0) = N, with flow time ϑ running from 0 to τ and the τ in the coefficients frozen. At fixed ε-order the flow converges 1/N-adically, so the x-adic and 1/N-adic expansions agree. Hence e^{βD} = 𝖲 + ε 𝖲_1 + O(ε^2), where 𝖲 f(N) = f(M), M := N + τ, and 𝖲_1 f = -ℓ f'(M), ℓ := log(M/N). Put ∇ := 𝖲 - 1. By (B2) and Step 3, modulo constants in y,

    𝖥 = x 𝖠 + 𝖡 + ε 𝖢 + ...,
    𝖠 = -N^3/(6σ) + κ_0 N,   κ_0 := (1+σ)/(24σ),
    𝖡 = -N^2/(4σ) - (1/24) log N + B_s,   B_s := Sum_{h>=2} σ^{h-1} Z_h N^{2-2h}.

Here ∂_τ is taken at fixed N, and primes are d/dN. Insert into (*).

Order ε^0: -2∇𝖠_τ - (1 + 𝖲)𝖠' + (∇𝖠')^2 = 0. This holds identically and is a check: the N^2, N^1 and N^0 coefficients each vanish, the last only because of the value of κ_0.

Order ε^1: 𝖢 does not enter, and

    -2∇(𝖡_τ + (σ'/σ)𝖠') - 2𝖲_1 𝖠_τ + (1 + 𝖲)𝖠'' + (1/N)(1 + 𝖲)𝖠' - (1 + 𝖲)𝖡' - (d/dN) 𝖲_1 𝖠 + 2(∇𝖠')[∇(𝖡' - 𝖠'/N) + 𝖲_1 𝖠'] = 0.

Use ∇𝖠' = -(M + N)/(2a_2). The pieces are:
- coefficient of ℓ: 2𝖠_τ'(M) + 𝖠''(M) - 2(∇𝖠') 𝖠''(M) = (MN - σ'/12)/σ^2;
- remaining source terms, not involving κ_0: (M + N)[σ'/(σ a_2) - 3/(2σ) - 1/(2a_2^2)] = -τ(M + N)/(2σ^2), from 2σ'τ - 3σ - τ^2 = -τ;
- terms with κ_0: (κ_0/a_2)(1/M + 1/N), using N/M - M/N = -τ(1/M + 1/N);
- with 𝓛[B] := -2∇B_τ - (1 + 𝖲)B' - ((M + N)/a_2)∇B', one has 𝓛[-N^2/(4σ)] = 0 identically, so the coefficient of N^2 is unconstrained by the next equation, and 𝓛[-(1/24) log N] = (1/(24 a_2))(1/M + 1/N).

Since κ_0 + 1/24 = (1 + 2σ)/(24σ) and σ/a_2 = τ, this gives

    𝓛[B_s] = (1/σ^2) { -ℓ (MN - σ'/12) + τ(M + N)/2 - (τ(1 + 2σ)/24)(1/M + 1/N) }.   (E1)

(E1) involves only the Z_h, no other unknown.

*Step 5 (solving (E1)).* [P] Let Bor: Q(τ)[[t]] -> N^{-1} Q(τ)[[1/N]], t^k -> k! N^{-k-1}, a bijection. Checking on monomials gives: Bor[e^{-τt} φ] = 𝖲 Bor[φ]; N Bor[φ] = φ(0) + Bor[φ_t]; (d/dN) Bor[φ] = -Bor[tφ]; Bor[(1 - e^{-τt})/t] = ℓ. Set 𝖯 := 1 - e^{-τt}, so that 𝖯_t = τ(1 - 𝖯) and 2𝖯_t + τ𝖯 = τ(2 - 𝖯). Write B_s = Bor[f_s/t], with f_s := Sum_{h>=2} σ^{h-1} Z_h t^{2h-2}/(2h-3)!. Then 𝓛[Bor[f/t]] = Bor[𝔏 f] for f in t^2 Q(τ)[[t]], where

    𝔏 f := 2𝖯 (f_τ/t - f_t/a_2) + (2 - 𝖯) f/a_2.

For the right side of (E1): with φ = 𝖯/t, MN ℓ - τ(M + N)/2 = Bor[φ_tt + τφ_t] = Bor[2𝖯/t^3 - τ(2 - 𝖯)/t^2], a regular series. Also 1/M + 1/N = Bor[2 - 𝖯].

Apply the formula 𝔏 to f_sing := -1/(σ t^2) + d_0, d_0 := (1 + 2σ)/(24σ). Using σ'/σ^2 - 2τ/σ^2 = 1/σ^2 for the coefficient of 𝖯/t^3, and ∂_τ d_0 = -σ'/(24σ^2), one gets

    𝔏 f_sing = (1/σ^2) [2𝖯/t^3 - τ(2 - 𝖯)/t^2 - (σ'/12) 𝖯/t + (τ(1 + 2σ)/24)(2 - 𝖯)],

which is exactly minus the Bor-preimage of the right side of (E1). Since Bor is injective, (E1) is equivalent to

    𝔏(f_sing + f_s) = 0.

On the other side, f_* := 1/(4 sinh(a_1 t/2) sinh(a_2 t/2)) = -1/(4 sinh(τt/2) sinh(a_2 t/2)) satisfies f_{*,τ}/(t f_*) - f_{*,t}/(a_2 f_*) = -coth(τt/2)/(2a_2). Therefore 𝔏 f_* = (f_*/a_2)[-𝖯 coth(τt/2) + 2 - 𝖯] = 0, because 𝖯 coth(τt/2) = 1 + e^{-τt} = 2 - 𝖯. Also f_* = f_sing + f_reg, with f_reg = Sum_{h>=2} β_{2h-2}(a_1, a_2) t^{2h-2}/(2h-3)!, since a_1 a_2 = -σ and a_1^2 + a_2^2 = 1 + 2σ. So W_0 := f_s - f_reg, in t^2 Q(τ)[[t]], satisfies 𝔏 W_0 = 0.

Uniqueness; **this is where the τ-symmetry is load-bearing.** (E1) alone does not determine the Z_h: the general solution of 𝔏 f = 0 is f_* g((1+τ)t) with g arbitrary, since the characteristics of ∂_τ/t - ∂_t/a_2 have invariant (1+τ)t. Suppose W_0 = w_n t^n + ... with w_n != 0 and n >= 2. The t^n coefficient of 𝔏 W_0 gives τ w_n' - (nτ - 1) w_n/a_2 = 0, so w_n = C τ^{-1} (1+τ)^{n+1}. But W_0 is invariant under τ -> -1-τ: σ and Omega are, since the weights τ and -1-τ swap, and f_reg is, since a_1 and a_2 swap. Invariance forces C[(1+τ)^{n+2} - (-1)^n τ^{n+2}] = 0, hence C = 0. So W_0 = 0, that is σ^{h-1} Z_h = β_{2h-2}(a_1, a_2) for every h >= 2. ∎

Sanity limits, not used. τ -> ∞ gives Int K Lambda^vee(1) = B_{2h}/(2h(2h-2)) = χ(M_h) (Harer-Zagier; consistent with PSM s.1.1(c)). τ -> 0 gives Int λ_h K = -(2h-3)! b_h, with b_h := Int λ_h ψ^{2h-2} > 0 (the λ_g value of Route A s.2.4 at j = 0).

Status. One referee pass (logic, literature and independent computation in one lens: sound with minor fixes, all applied above). The referee confirmed the value by an MV-free computation (Mumford/Faber recursion and DVV) at h <= 4 at three τ and at h = 5 at one τ [E]. Scope: Step 1 settles (B1) for all g and μ, the same derivation as Proposition MV(A) of the MVmatch report; it does not by itself establish Lemma A_P2(b), whose status is given in s.7.4.


## 8. Status of every statement; remaining inputs

### 8.1 Status table

    Statement                                        Status      Rests on / source
    ------------------------------------------------ ----------- ------------------------------------------------
    Lemma 1.2(a),(b) (t_c, ρ, ρ^2 Delta_0)           P           algebra
    Lemma 1.2(c) (t_c = t_cf)                        P mod C     BS 5145, eq:pfinf (read by round-4 critics);
                                                                 used only for Corollary G
    Lemma 1.4 (code bookkeeping)                     P           reading hv_lib3.m, hv_lib4.m, lp2_lib.m, lp2_percolour.wls
    Theorem A, Cor 2.2 (PF, uniqueness)              P           algebra; R3 Part I; Theorem A re-derived by round-5 referee
    SP1-SP5 (formal stationary phase)                S (written) R3 Part C; interface only
    Lemma 2.1 (critical points)                      P           R2:Lemma 1; interface only
    Lemma 2.3 (analyticity)                          P           R2:Lemma 2; interface only
    Lemma 2.4 (PF under the integral)                P           SP2, SP3; interface only
    Lemma 2.5 (WKB uniqueness)                       P           R2:Lemma 4; interface only
    Prop 2.6 (normalisation; eq:normNLR derived)     P mod C     (C1), (C2); R2:Prop 5; interface only
    Thm 2.7 (continuity, closed form, Rec = Wick)    P mod C     (C1)-(C3); R2:Thm 6, Prop 8; R3:Cor B; interface only;
                                                                 reading of BS 5257-5267 re-checked by referee
    Cor 2.8 = Lemma CL (corrected)                   P           algebra; rationality at colour 0, frame B only
    Cor 2.9: F2, L1, Form, R, Z-infinity             P           closed form + Theorem A
    Prop 2.10 (symplecticity, K = 1, orthogonality)  P           Theorem A only; NEW in round 5, one referee pass
    Prop 2.11 (Givental-Teleman on T_F)              P mod C     (C1)-(C3), (C5); not used for Theorems LC_P2, R_P2
    Prop 2.12 (lem:modrel, both relations)           P           no citation; NEW in round 5, one referee pass;
                                                                 interface only
    Lemma D0 (L1, F3 for Dil_0)                      P           Cor 2.9; written here
    Prop 3.1 (closure, general c1)                   P           Theorem A; not used (its c1 = 1 form is S1 of Prop 2.10)
    Lemmas 4.1-4.3, Prop H2, Lemma G                 P           s.3.1, Prop 2.10
    Lemma R, Theorem E00                             P           s.3.1, H2 Steps 1-2, Prop 2.10
    Theorem M                                        P           Cor 2.9(R), Prop 2.10, H2
    Theorem LC_P2                                    P           for the decoration sum as defined; no BS input
    Fpot_h on 𝔇, (H0)                                P           dilaton equation; NEW in this edit, unrefereed
    (S_P2), (C_P2), (L0_P2)                          P           string/dilaton equations
    Lemma 7.1 (range)                                P           stability
    Lemma A_P2(a),(c)                                P           degree count
    Lemma A_P2(b) (vanishing)                        C           lemmaA/LEMMA_A_PROOF.md Theorem L, P there modulo (MV) and
                                                                 textbook facts; becomes P mod (MV) after its open point O1;
                                                                 see s.7.4 and 8.2
    Theorem R_P2                                     P mod C     Theorem LC_P2, Lemma A_P2(b) for h <= g-1
    Theorem VI (vertex identity (VI))                P mod C     LLZ eqs. (1), (5), boson-fermion, string/dilaton, (7.2);
                                                                 NEW in round 5, one referee pass
    Corollary SV_P2                                  P mod C     above + Lemma A_P2(b) at h = g, Theorem VI
    Proposition 13.2_P2 (G_g = GW^cf_g)              O           s.8.3
    Corollary G (BS rmk:cfgap)                       Conditional on Proposition 13.2_P2. NOT CLAIMED.

"Interface only" means: used for Proposition 2.11, Proposition 13.2_P2 or Corollary G, and not in the proofs of Theorems LC_P2, R_P2 or Corollary SV_P2. Round-4 critics asked for H2 Step 3, E00 and (M-d) to be relabelled P mod (C1)-(C3); since Proposition 2.10 is now P, that relabelling is unnecessary. All of these theorems are about the decoration sum (0.1) as defined.

### 8.2 Remaining inputs for Theorem R_P2 (and Corollary SV_P2)

The following is the outcome of a dependency reading, not of a full referee pass over the document.
1. **Lemma A: LEMMA_A_PROOF.md (s.7.4).** Literature: (MV), the Marino-Vafa formula, LLZ eq. (1) (arXiv v3 read, JDG numbering not confirmed), with OP Theorem 1 as a second published proof; and the textbook facts listed in s.7.4. The 2D Toda/mKP equations are not a citation in the proof of record, since the P_k are derived in s.4 of that file. The matching of the MV series with Xi_h is no longer an open item (that file's Proposition 2.1, s.3 and Proposition 8.1). Open there: O1, the assembled text has not had its own referee; O2, citation details, none mathematical; O3, the proof of (MV) in LLZ was not re-checked.
2. **Theorem VI, for Corollary SV_P2 only.** LLZ eq. (1) and Theorem 1 eq. (5); the boson-fermion correspondence (vertex-operator form of ψ, ψ^*); the principal specialisation. The string and dilaton equations and (7.2) are standard or proved here. The ABDKS/Eynard route is no longer used.
3. **BS.** No BS input enters Theorems LC_P2 and R_P2, or Corollary SV_P2 on the recorded route. This rests on the dependency reading of s.3-s.7 by the round-5 modrel author and referee; the η^{ab} and the operator eq:pf enter as definitions. Identification with GW still needs (C1)-(C5) and Proposition 13.2_P2, s.8.3.

Unconfirmed citations are collected in s.7.4 and s.7.8: JDG numbering of LLZ; theorem numbers for the mKP bilinear identity in Date-Jimbo-Kashiwara-Miwa and Miwa-Jimbo-Date; Macdonald I.3 Ex. 2. F2 has **not** been taken from thm:cfasym: since round 2 it follows from the closed form, and eq:Rcf is then a consistency check.

### 8.3 Frame-B formula on T_F (Proposition 13.2_P2) and the gap corollary as the final step

**Proposition 13.2_P2 [O].** For g >= 2, Cf(GW_g^{T_F}) equals G_g of (0.1).

What is available:
- Proposition 2.11 (Givental-Teleman on T_F, mod (C1)-(C3), (C5));
- Theorem 2.7(e) (R^{T_F}_{0i} = h_w × frame-B row, with z-only factor exactly 1);
- Lemma 1.2(c) (t_c = t_cf, A_cf = C1cf);
- Proposition 2.12 (lem:modrel) and its dictionary, identifying Definition 1.3 at c1 = A with eq:RHb;
- the Cf convention, s.1.3;
- N1, N2, N5 of the round-1 frame-B skeleton are resolved by Proposition 2.6, Theorem 2.7, Definition 1.3.

What remains (from the round-1 skeleton, which has **no verifier**; the pathway note warns that it is a theorem, not a paragraph):
1. *Frame B.* Omega_w = h_w·1 by Mumford, and the PPZ group property gives R' = R h_w^{-1}: port of conifold Prop 13.2(e). In particular one must show that the code's dilaton leaf q_t = -ψ S (s.1.4) is the Givental translation of R', with the scalar translation of h_w absorbed by Omega_w. The scoping note warns that on T_F the mechanism of conifold Lemma 13.3 (the h_nw cancellation) has no direct analogue, since the large-radius point of T_F is not semisimple. The route through T_max and continuation (Proposition 2.11) has to be written out with the frame-B twist.
2. *N4: closed.* lem:modrel, including the second relation E(y) = (D - X)/2, is proved in Proposition 2.12, and it supplies the dictionary Δ_c A^2 = 1/D, d2 = e_1(α*) E(y) identifying Definition 1.3 at c1 = A with eq:RHb.
3. *BS line ~4056.* Regularity of GW_g^{T_max} on the diagonal, argued only briefly (properness); used in Proposition 2.11. To be checked.
4. *N3.* The lift of identities to G_{T_F}[A^{±1}, F], using algebraic independence (C4). Allowed, needs writing. Since Cf is defined as a ring homomorphism, it acts graph by graph with no on-shell input.
5. *Verification.* An independent referee pass on the whole frame-B skeleton.

**Remark on N4 and Proposition 2.10.** In round 4, Proposition 2.10 step 3 compared with BS's rows at c1 = A, which made lem:modrel a load-bearing input of Theorem LC_P2, and the draft wrongly attributed c2 = D to eq:ADEFXdef (it is the first relation of lem:modrel). Both issues are gone: Proposition 2.10 is now proved from Theorem A, and lem:modrel is Proposition 2.12. The audit of BS's own use of lem:modrel is Remark 2.13.

Round-1 items N1-N8 are recorded in journal wf_e928f771-9e3 (only N4 appears in ROUND1_STATUS_2026-11-02.md).

**Final step (gap corollary).** Once Proposition 13.2_P2 is proved:
- Theorem R_P2 + Corollary SV_P2 + Lemma 1.2(c) give GW^cf_g = (-3)^{g-1} β_{2g-2}/t_cf^{2g-2} + O(1), i.e. BS rmk:cfgap, the full refined conifold gap.
- No TR for the target theory enters. Watch-points (round-1 skeleton): never use BS rmk:refmm (Eynard spectral curve), BKMP/FLZ remodelling for K_P2, or HK10 / Coates-Iritani 2018 / BFGW Conj 1.9 as input, as they state or anticipate the gap.
- A separate referee pass is required.


## 9. Errata and imprecisions found in BS v2 (BS_v2_2026-11-02/main.tex)

None is load-bearing for the proofs above, after the fixes indicated.

Errors:
1. **Line ~3382-3384.** W_i on T_F is printed as -(E/9)(3Y+E)^2. It should be +(E/9)(3Y+E)^2, by eq:Wadm (3225), eq:delta1 (2971) and direct computation. The code uses the + sign: sqDsym^2 = ε1 ε2 E(E+3L)^2/3 = 3 ε1 ε2 W. (The unverified round-1 skeleton quotes the wrong sign; do not copy.)
2. **eq:crituiWi (4391).** In q^{1/3}-coordinates the determinant is 243·W, i.e. off by a constant factor in that chart.
3. **Line ~4370.** "q_0^cr = c(1, ω, ω^2)" cannot be a critical point in q-coordinates: q_i ∂_i 𝒲 = q_i - α + P/3 forces all q_i equal. It is a vector of q^{1/3}-coordinates, mapping to q = c^3(1,1,1).

Imprecisions and gaps, with the fix here:
4. **BS 4337-4343** present the y -> 0 limit of Asym_i "in particular" after eq:oscintP + eq:normNLR. It cannot then be used to derive eq:normNLR. Fix: the direct decoupling of Proposition 2.6.
5. **BS ~4344-4346** assert regularity at α_i = E/3 and eq:oscintP "by continuity". Fix: Theorem 2.7(a).
6. **BS 3497-3510 (Normalisation)** is asserted with five-weight constants without proof. For T_F purposes, fixed by Proposition 2.6 modulo the BS definition of the canonical R-matrix.
7. **BS, proof of thm:cfasym (5301-5334).** The s_3 and s_3^2 insertion identities are stated only at leading order in ξ, and the row-2 insertion (~5320) is frame-dependent: Cf(R_2) involves F_cf, E_cf. They must not be read as exact. Round-1's "Route B" is dropped. Not needed here.
8. **sec:cfgen (5140).** Cf is described through τ -> -1/(3τ), which does not fix the Hauptmodul y. Cf must be defined as the identity on G_{T_F} (BS's commented-out identification Y^cf = Y) and extended to A^{-1}. Here this is a stated convention (s.1.3).
9. **prop:fingenR (3611)** is stated for rational tori only, but is used at T_F. Fix: Theorem 2.7(e).
10. **exmpl:stirling (4260)** is analytic, for α in R^+, and is used at α = -E. Fix: formal SP5.
11. **lem:modrel (2695)** E = (D - X)/2: BS's proof is a sketch ("We omit the details"), and the first relation is cited. Fix: Proposition 2.12, whose proof can replace 2697-2749 verbatim. Also: eq:SHb (about 3299) cites eq:allc, whose printed row 2 already has modrel applied; it should cite eq:ssH1 (Remark 2.13(a)). BS's derivation of Delta_i^{-1} passes through modrel-simplified formulas; a modrel-free derivation is in Remark 2.13(c).
12. **Line ~4056:** the regularity of GW^{T_max} on the diagonal is argued briefly.
13. **Wick numbers.** Only r_0, r_1, r_2 are displayed (5257-5262); the recurrence (Rec) is not in BS. It is proved here (Theorem 2.7(e), step 3).
14. **Line references.** eq:pf label at 2628 (displayed about 2627-2632); eq:Wadm at 3225; eq:Rcf at 5334; thm:cfasym at 5301; lem:cfpole at 5407; rmk:cfgap at 5422; eq:Ifun at 2621; eq:sympl at 3084; eq:SHb at 3305; eq:RHb at 3573; prop:ResConifold at 1708 ("Prop. 3.1" is the numbering of arXiv:2410.00118, not a v2 label); prop:cfgraphsum at 6366. Earlier round notes quote slightly different numbers.


## Appendix. Computer checks (evidence only; none is a proof step)

The round-4 critic's objection, that the items S1, S3 and proof (ii) of S2 made finite-order evidence load-bearing, no longer applies: those steps now have written all-order proofs (s.7.4, s.7.8), and gap_identity_check.pl, dexrel_check.pl and the data comparisons are demoted to evidence.

Scratchpad root: SP := <user home>/AppData/Local/Temp/claude/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/scratchpad. Project root: LP := G:/My Drive/Projects/Old/refGW/hodge_vertex/lp2_explore. Copies of the round 1-3 scripts are in LP/round1_scripts, LP/round2_scripts, LP/round3_scripts. Scratchpad paths may be purged.

**Before round 1 (LP; RESULTS_2026-11-01.md).** Exact Perl/Math::BigRat and Mathematica (main session).
- hrow.pl, hrow_out_*.txt: H-row form through t_c^7, z^9 at 4 weight pairs.
- edges_sing.pl, edges_sing_out_*.txt: nine edges regular; rank one at t_c = 0; loop and leg support, m+n <= 8, s^4; XX(0) = -2/9 needed.
- lp2_percolour.wls, lp2_percolour_g2_out.txt, lp2_percolour_g3_out.txt: genus 2 (6 edge graphs, 2 weights) and genus 3 (41 edge graphs, 2 weights), all "REGULAR colouring by colouring"; in the single-vertex graph only the singular colouring has a pole, which equals the gap.
- lp2_vertex_conjecture_g2.wls, _g3.wls: leg and loop configurations have no polar term; 42 heavier genus-3 configurations were skipped (they need rows beyond z^9). Readout caveat in RESULTS item 4.
- frob_exact_out.txt, frob_vs_row_out_*.txt: Frobenius cross-check; only simple poles at z = -E/j.

**Round 1 (SP/lp2_round1/, folders cfasym, cfasym_verify, closure, closure_verify, diff, frameB, h2, h2_verify, regular, regular_verify, residue, residue_verify).**
- row_vs_closed.pl: stored cPP_z9 equals the closed form at finite rational L, z^2..z^9, 4 weight pairs, z -> -z control fails.
- closure: A_1..A_9 = 0 to s^4; off-shell anomaly A_2 = -c1 Rd E_21.
- h2: no ψ-linear term in σ_2 through t_c^7.
- residue: Lemma R reproduces the code's E_00 exactly, 45 entries, s^-2..s^6, 3 weight pairs; conifold control 644/0; off-shell controls fail.
- dexrel_check.pl: E(y) = (D - X)/2 to y^40. Evidence only since Proposition 2.12.

**Round 2 (SP/lp2_round2/: H1_bs, H1_code, continuity, referee_code).**
- pfcheck.pl: stored rows satisfy eq:pf exactly through z^10 at 6 points (cannot detect a wrong z-only factor by construction).
- closed form satisfies q^3 𝒪 f = 0 exactly through z^19.

**Round 3 (SP/lp2_round3/: h2, H2_ref2, E0j, E0j_ref2, PF, PF_ref1, PF_ref2).**
- H2 and M statements at finite orders; consistent off-shell section c1 = -wc(1 + Q Δ_c) at (7,3), Q = 1/10 and (4,9), Q = -3/7 (k = -1..5, values indexed from k = 0).
- Direct rank-one check of Ehat_0j at t_c = 0 in all entries.
- Recurrence for r_l against Wick numbers; PF operator identity.

**Round 4 (SP/lp2_round4/).**
- thmR/checks.py: (S_P2), (L0_P2) on the genus-0 trivial CohFT at E = 7/3 (control of the algebra); control J/E instead of J/E^2 fails. Also the single-vertex pole against the logged lp2 values: 5 of 5 cases (g = 2 at x = 1/3: -83/6480, x = 1/8: -9773/327680; g = 3 at x = 1/3: -169/27216, x = 1/8: -1666067/117440512, x = 4/5: -1414951/224000000).
- thmR_ref2/ref2.py: (S_P2), (L0_P2) at genus 0 and genus 1, referee's own recursion from ⟨τ_1⟩_1 = 1/24, 11 cases; the controls (J/E; D instead of D+2 at h = 1) fail.
- thmR_ref2/zinf.py: L2b(θ S_∞) = 0 through V^14; perturbed r_2 control non-zero.
- lemmaA/gapcheck.py: (-3)^{g-1} gap_g against stored local-P2 values at g = 2, 3, 4: True. The stored values (hodge_vertex/README.txt lines 28, 70, 99; hodge_vertex_genus3_summary.txt; hodge_collect_genus4.txt) are total leading polar coefficients. They agree with the single-vertex pole at leading order, since edge graphs are subleading; the genus-4 per-graph line for graph {4} shows the gap followed by zeros.

**Lemma A data (hodge_vertex_conifold/).** hvc_multiprobe_out.txt, hvc_multiprobe_g4_out.txt, hvc_multiprobe_g4b.wls (M_irr(3), M_irr(4)); lemmaA/E_formulas/mirr_g7.txt (all 194 genus-7 points zero, 77 controls non-zero); lemmaA/round2/critic/gap_identity_check.pl (44 comparisons through h = 12).

**Round 5 (SP/lp2_round5/). Exact rational arithmetic, or F_p where stated. Checks by the MVmatch and ThmL owners and referees are recorded in LEMMA_A_PROOF.md and journal wf_9d7fd94e-8d7.**
- modrel/modrel_series.py: L2b A = 0, (K), L2b P_I = 0, A^2 D = X, 2E(y) = D - X and eq:limitconv through y^60; eq:Ian coefficients reproduced from eq:Ifun for d <= 6 at a generic rational α.
- modrel/pairing_check.py: Π_ij = δ_ij on the closed form alone, all i <= j, through z^6 at (1,1), (5,2), and through z^5 at (7,-3), (1/2,3/5). Controls (r_2 perturbed, η^00 doubled) fail.
- modrel/closure_frob_check.py: d2 = e_1 E(y) and Δ_c c1^2 = 1/D at c1 = A through t^24 (Laurent truncation); Definition 1.3 = eq:RHb through z^4 up to t^4 at colour 0 and (5,2), (1,1); floating-point test of the modrel-free Phi_2i.
- modrel_ref1/thm3_fp.py (F_p, p = 2^61 - 1, random section c1, random t0): Π_ij = δ_ij for all colour pairs through z^14 at 9 weight pairs; Proposition 3.1 for general c1 through z^13. Controls (d2 doubled, r_3 perturbed, η^00 doubled) fail.
- modrel_ref1/thm1_direct.py: I-function expanded directly from eq:Ifun at 4 α, including T_F; A^2 D = X, 2E(y) = D - X, eq:DEXrel through y^22; (K) through y^200. Control fails.
- modrel_ref1/generic_alpha.py (F_p): at generic α, η_ab from Remark 2.13(c), and an independent WKB solution of eq:pf with rows from eq:RHb: B_ij = 0 for i != j, θ B_ii = 0 through z^5. Controls fail.
- VI/vi_check.py: (*) on the Schur side, p_1^n, n <= 7 at τ = 2/7, -5/3; Z_h extracted from MV data, h = 2..6, equal to β_{2h-2}; (E1) with dual numbers in τ.
- VI/vi_grade.py: (*) in monomials x^a w^b, x-order <= 12: grade <= 1 vanishes with Z_h = β_{2h-2}/σ^{h-1}. Perturbing Z_3 gives non-zero grade 1.
- VI_ref1/faber.py, h5.py: MV-free values (Mumford/Faber recursion and DVV) of σ^{h-1} Z_h = β_{2h-2}(-τ, 1+τ) for h = 2, 3, 4 at τ = 3/11, -7/4, 5/2, and for h = 5 at τ = 3/11. Control β·σ fails.
- VI_ref1/lastchecks.py: slice normalisation (B1) against Faber integrals for (h, a) up to (3, 2); (E1) with no mismatch for N^1 down to N^{-33} at τ = 3/11, -7/4, 13/5. This closes the top-order flag of VI/vi_check.py as a truncation artefact. Kernel control: adding (5/3) f_* ((1+τ)t)^4, which breaks τ-symmetry, still solves (E1), so the symmetry step is genuinely needed.
- VI_ref1/schur_modp.py (mod 2^61 - 1): (*) on p_1^n, n <= 12; Z_h for h = 2..8 at τ = 3/11, -7/4, 13/5 by the w-polynomial fit, equal to β.
- VI_ref1/grade.py: with B_s and higher grades as arbitrary random data, grade 0 of (*) vanishes and grade 1 equals 𝓛[B_s] - RHS(E1) exactly. Controls: perturbing the log coefficient or the genus-1 constant fails, and perturbing the -1/4 of G_0 changes nothing, as 𝓛[N^2/σ] = 0 predicts.
- draft/check_transcription.py: this text's formulas for (E1), 𝔏 and the (H0) base-point identity, re-evaluated exactly with the notation of s.7.8.
