# Round-5 patch, part A: header, s.0, s.1. Entries: (start, end, prefix of original line `start`, new text); 1-indexed inclusive.
PA = []
PA.append((1, 3, '<!-- Round 4', r'''<!-- Round 5 (2026-11-02): edited in place from the round-4 consolidated draft (backup: LP2_REGULARITY_PROOF_draft_round4.md). Applies ROUND4_CRITICS_2026-11-02.md item by item, integrates the round-5 reports "modrel" and "VI" with their referee fixes, and points Lemma A to hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md. Built on ROUND2_PROOFS_rows_on_TF.md, ROUND3_PROOFS_H2_edges_PF.md and the status notes. Journals (under <user home>/.claude/projects/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/subagents/workflows/<run>/journal.jsonl): round 1 wf_e928f771-9e3 (this is the only record of the round-1 frame-B skeleton and of the items N1-N8 quoted below, except N4), round 2 wf_1ab43e48-d73, round 3 wf_74ed0422-50f, round 4 wf_bab37015-382 (reports thmR, lemmaA and their referees), round 5 wf_9d7fd94e-8d7. NOT yet refereed as a whole. -->'''))

PA.append((10, 10, '- **[P]**', r'''- **[P]** proved in this text, to all orders. Where a proof is taken over from an earlier round, the round and item are given (R2:x = ROUND2_PROOFS_rows_on_TF.md item x; R3:x = ROUND3_PROOFS_H2_edges_PF.md item x; the round-4 reports thmR, lemmaA and their referees are in ROUND4_CRITICS_2026-11-02.md and journal wf_bab37015-382; the round-5 reports modrel, VI, ThmL, MVmatch are in journal wf_9d7fd94e-8d7).'''))
PA.append((12, 12, '- **[S]**', r'''- **[S]** standard and written out here with its proof: formal stationary phase (SP1-SP5), s.2.2. An [S] item is not an open input.'''))

PA.append((16, 16, '**Rules', r'''**Rules (user rulings, ALLGENUS_PATHWAY_2026-11-01.md s.6).** Allowed: BS genus-0 mirror and oscillating-integral statements, cfasym, finite generation and algebraic independence; KP/Toda; topological recursion (TR) on the Hodge side. Forbidden: TR for the target theory, including the local-P2 mirror curve. Everything must hold to all orders: finite-order computations are evidence only and may not fix a formula or normalisation that a proof relies on. The gap corollary is in scope as a final step. No step below uses TR for the target theory. Hodge-side TR enters only through the backup route for Lemma A (s.7.4) and the cross-check of Theorem VI (s.7.8), and neither is a proof step.'''))

PA.append((18, 18, 'Contents.', r'''Contents. 0 Summary. 1 Conventions. 2 Rows on T_F, including s.2.7 (lem:modrel and the interface with BS). 3 The singular colour: F1, F2, L1, Z-infinity, closure, K = 1, Dil_0. 4 Proposition H2. 5 The residue lemma and E_00. 6 Regular colours and mixed edges (Theorem M). 7 Hodge side and Theorem R_P2, including s.7.8 (proof of the vertex identity). 8 Status of every statement; remaining inputs. 9 Errata in BS v2. Appendix: computer checks.'''))

PA.append((25, 38, 'Fix ε1', r'''Fix ε1, ε2 with ε1 ε2 E != 0, where E := ε1 + ε2. In the code these are called e1, e2. The Hodge weights are w = (E, -ε1, -ε2). For g >= 2, the frame-B Hodge-vertex graph sum of genus g is

    G_g := Sum_Gamma |Aut Gamma|^-1 Sum_{c: V -> {0,1,2}} T(Gamma, c),
    T(Gamma, c) := Sum_dec prod_v [ Delta_{c_v}^{g_v - 1} (1/m_v!) Int_{Mbar_{g_v, n_v + m_v}} Omega_w prod_{h at v} ψ_h^{a_h} prod_{k=1}^{m_v} ψ_{n_v+k}^{d_{v,k}} · prod_k Dil^_{c_v}(d_{v,k}) ] prod_e Ehat_{c c'}(a_e, b_e).      (0.1)

The ingredients are as follows.
- Gamma runs over connected stable graphs of genus g without legs; n_v is the number of half-edges at v.
- **Decorations; (0.1) is their definition.** A decoration assigns an exponent a_h >= 0 to each half-edge (a_e, b_e denote the exponents at the two half-edges of e), and to each vertex an integer m_v >= 0 together with an ordered tuple (d_{v,1}, ..., d_{v,m_v}) of dilaton degrees >= 2, with weight 1/m_v!. Since Omega_w has mixed degree, the integral vanishes unless Sum a_h + Sum d_{v,k} <= dim = 3g_v - 3 + n_v + m_v. As d >= 2 this forces m_v <= 3g_v - 3 + n_v. So the decoration sum is finite, and T(Gamma, c) is a finite sum of products of the factors below; this is what gives "analytic at the conifold" its meaning in Theorem R_P2.
- Omega_w is the normalised Hodge class prod_{a=1}^{3} Lambda^vee_{g_v}(w_a)/w_a^{g_v} (s.1.5). It is the same at every colour and independent of y.
- Delta_i, Dil^_i and Ehat_ij are the normalised frame-B factors of s.1.4. They are built from the rows of Definition 1.3, taken in the conifold polarisation c1 = C1cf.

Colour 0 is the singular colour: on the real branch t = y^{1/3} -> -1/3 it is the unique colour whose critical point escapes (Lemma 2.1).

(0.1) is the definition of the decoration sum. It reproduces the code's value (Lemma 1.4). Theorems LC_P2 and R_P2 are statements about this sum. Identifying G_g with a Gromov-Witten potential is a separate matter (Proposition 13.2_P2, s.8.3).'''))

PA.append((42, 42, '**Theorem LC_P2', r'''**Theorem LC_P2 (leaf structure at the singular colour; F1-F3, L1-L3, to all orders in t_c and z).** [P: no BS input, since Proposition 2.10 is now proved from the closed form]'''))

PA.append((55, 55, '**Theorem R_P2', r'''**Theorem R_P2 (regularity of graphs with an edge).** [P mod Lemma A(h), 2 <= h <= g-1]'''))

PA.append((59, 59, '**Corollary SV_P2', r'''**Corollary SV_P2 (the single-vertex graph).** [P mod Lemma A(h) for 2 <= h <= g, and Theorem VI, which is P mod LLZ eqs. (1), (5) and the boson-fermion correspondence (s.7.8)]'''))

PA.append((65, 65, 'By the Hodge-side', r'''By the vertex identity (VI), Theorem VI of s.7.6, this equals (-3)^{g-1} β_{2g-2}(ε1, ε2) t_c^{2-2g}, where β_n := (n-1)! [t^n] 1/(4 sinh(ε1 t/2) sinh(ε2 t/2)) (BS eq. at line 959). Hence the polar part of G_g is exactly (-3)^{g-1} β_{2g-2} t_c^{2-2g}.'''))

PA.append((72, 80, '- **BS genus-0 inputs', r'''- **BS genus-0 inputs (allowed).** Since round 5 these are **not** inputs of Theorems LC_P2 and R_P2, or of Corollary SV_P2 on the recorded route for (VI): Proposition 2.10 is proved from Theorem A alone, and s.3-s.6 contain no other BS link (dependency reading by the round-5 modrel author and referee). (η^{ab} and L_pf enter those theorems as definitions.) The BS inputs are needed only to identify the rows and the sum with BS's objects, i.e. for Lemma 1.2(c), Proposition 2.6, Theorem 2.7, Proposition 2.11, Corollary G and Proposition 13.2_P2. They are labelled (C1)-(C3) as in R2:s.0; the parenthesised labels (Cn) are never to be read as C1cf or C1i.
  - (C1) eq:Ifun 2621, eq:J=I 2663, eq:asymp 2943 (with about 2950) and eq:canR 3566: for generic α, e^{u_i/z} R^α_{0i} formally solves eq:pf (2628).
  - (C2) (Normalisation) 3497, eq:phi 3031 and the weights at about 2383: the value at y = 0 and Psi_{0i} = Delta_i^{-1/2}.
  - (C3) eq:RHb 3573, eq:ADEFXdef 2681, eq:etaflat 2419 and (Symplecticity) 3496. lem:modrel 2695 is no longer cited: both of its relations are proved in Proposition 2.12.
  - Also used: (C4) prop:algindABC 2854, only to make Cf well defined; and t_cf at line 5145 with eq:pfinf 2757, only for Lemma 1.2(c).
  - For Proposition 2.11 and Proposition 13.2_P2 only: (C5) eq:gtH 3799, thm:haegen 3902, lem:lho 3247 and the regularity on the diagonal at about line 4056.
- **Lemma A (Hodge side).** This is the conifold Lemma A, transferred in s.7.4. Its written proof is hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md: Theorem L for every h >= 2 and every non-empty R, proved there modulo the Marino-Vafa formula (LLZ eq. (1); OP Theorem 1 is a second published proof) and textbook facts. The assembled file has not yet had its own referee; see s.7.4 and s.8.2.
- **Hodge-side standard facts.** These are the string and dilaton equations for Omega_w; they hold because Lambda^vee is pulled back under forgetful maps.
- **Vertex identity (VI).** Used only for Corollary SV_P2. It is now proved in all genera (Theorem VI, s.7.8) from LLZ eqs. (1) and (5), the boson-fermion correspondence, the string and dilaton equations and (7.2). It no longer rests on the ABDKS/Eynard normalisation.'''))

PA.append((86, 86, '- A literature-free', r'''- A literature-free proof of Lemma A or of Theorem VI: both are proved modulo the literature listed in s.8.2.
- That the whole text has been refereed. Proposition 2.10, Proposition 2.12 and Theorem VI have had one referee pass each (sound with minor fixes, applied here); the extension of Fpot_h in s.1.5 and (H0) of s.7.1 are new in this edit and unrefereed.'''))

# s.1.1
PA.append((100, 100, '  Round-3 referee caution', r'''  Caution: the script edges_sing.pl (line 129, variable $vi) uses v for 1 + 3t_i. That is a different quantity, and it is never used here.'''))
PA.append((102, 104, '- **Edge variable.**', r'''- **Row and edge variables.** z is the argument of a row. ψ, ψ', w are edge variables: a half-edge with variable ψ carries the row at z = -ψ, i.e. P^^{(i)}(-ψ). So Ehat_ij(ψ, w) of (1.1) has edge variables as arguments, not row arguments.
- **Leaves.** λ_j(ψ) := 1/(E - jψ) for j in Z \ {0}, and λ_0 := 1/E. Further Λ_m := span{λ_1, ..., λ_m}, Λ_0 := 0, Λ_m^+ := C + Λ_m; in bounds, Λ_m is read as Λ_{max(m,0)}.
  **Normalisation warning.** R2:Cor 9 uses (1 - jψ/E)^{-1} = E λ_j, while R3 and this text use λ_j. Spans are unaffected, but the coefficient statements σ_0^(0) = E λ_1 and Sum_j β_j/j = E δ_{k0} are correct only in the present normalisation.'''))
PA.append((106, 107, '  - (LAM1)', r'''  - (LAM1) -ψ λ_j = (1 - E λ_j)/j, for j != 0.
  - (LAM2) (E - Nψ) λ_j = N/j + (1 - N/j) E λ_j, for j != 0; in particular (E - Nψ) λ_N = 1 for N != 0.'''))

# s.1.2 remark
PA.append((125, 125, 'Remark: (c) is new', r'''Remark: (c) was new in round 4. Critics 1 and 2 of round 4 read eq:pfinf (it is L2b∘θ) and BS 5145 and agree with the proof, but it has had no separate referee pass. It is used only for Corollary G. The normalisation of t_c enters every exact constant through ρ(0) = 1, and through Z-infinity (Corollary 2.9, s.3.1). All filtration statements hold for any ρ(0) != 0.'''))

# s.1.3
PA.append((129, 129, '- Θ(x) :=', r'''- Θ(x) := exp Sum_{m>=1} B_{2m} x^{2m-1}/(2m(2m-1)), so Θ(-x) = 1/Θ(x). The two-variable version used in SP5 and Theorem 2.7 is BS's (line 4278): Θ(z, x) := exp[-Sum_{k>=1} B_{k+1}(x) (-z)^k/(k(k+1))], with B_n(x) the Bernoulli polynomials. Then Θ(z, 0) = Θ(z), Θ(z, x+1) = (1 + xz) Θ(z, x), Θ(z, -m) = Θ(z)/prod_{j<=m} (1 - jz);'''))
PA.append((134, 134, '- S^[m](v, z)', r'''- S^[m](v, z) := Sum_l r_l (vz)^l / prod_{j=1}^{l+1-m} (1 + jz/E), where a product with upper limit < 1 is read as 1 (the case l = 0, m = 2 occurs in Lemma 4.1).'''))
PA.append((147, 149, '- With c1 = A', r'''- With c1 = A these are eq:RHb (3573), through Δ_c A^2 = 1/D and d2 = e_1(α*) E(y). Both are lem:modrel, proved in Proposition 2.12 (the first is its first relation, not a definition). With c1 = C1cf = A_cf, and Cf taken **by convention** to be the identity on G_{T_F} (well defined by (C4); BS sec:cfgen describes Cf only through τ -> -1/(3τ), which does not fix the Hauptmodul), they are Cf(R_{ai}). A general solution of L2b gives the c1-polarisation, not Cf. BS's A satisfies L2b because A = θ of a period of L2b∘θ.

Theorem 2.7 identifies Definition 1.3 with the T_F limit of BS's canonical rows (R2:Definition 7); the z-only factor is exactly 1. Since round 5 that identification is **not** used in the proofs of Theorems LC_P2 and R_P2 or of Corollary SV_P2: everything in s.2.6 and s.3-s.6 is algebra on the closed form, including symplecticity (Proposition 2.10). The identification is needed for Proposition 2.11 and Proposition 13.2_P2.'''))
PA.append((151, 151, '**Pairing.**', r'''**Pairing.** ⟨a, b⟩ := Sum η^{ab} a_a b_b, with η^00 = ε1 ε2 E^3/9, η^11 = -3 ε1 ε2 E, η^12 = η^21 = -3 ε1 ε2, all others 0. Note η^11 = E η^12 and η^22 = 0. In Theorems LC_P2 and R_P2 this is a definition. It agrees with lp2_lib etainvF and, for the last three values, with BS line 5395; matching it with eq:etaflat is needed only for the interface (s.2.7).'''))

# s.1.4
PA.append((155, 155, '- **Dilaton factor:**', r'''- **Dilaton factor:** Dil^_i(d) := (-1)^d [z^{d-1}] P^_0^{(i)} for d >= 1, equivalently q_t^{(i)}(ψ) := Sum_{d>=1} Dil^_i(d) ψ^d = -ψ S(v_i, -ψ). Since [z^0]S = 1, we have Dil^_i(1) = -1 at every t: the d = 1 term is the dilaton shift, and the decorations of (0.1) have d >= 2. **Dil^_i uses only component 0 of the row.** So it does not depend on c1, and neither Lemma G nor the sign of c1 is involved. We write q_t := q_t^{(0)}.'''))
PA.append((156, 160, '- **Edge factor:**', r'''- **Edge factor:** with edge variables ψ, w (rows at z = -ψ, z = -w; s.1.1) and P_ij(ψ, w) := -⟨P^^{(i)}(-ψ), P^^{(j)}(-w)⟩,

      Ehat_ij(ψ, w) := [P_ij(ψ, w) - P_ij(-w, w)]/(ψ + w),   Ehat_ij(m, n) := [ψ^m w^n] Ehat_ij.        (1.1)

  By symplecticity (Proposition 2.10), P_ij(-w, w) = -⟨P^^{(i)}(w), P^^{(j)}(-w)⟩ = -δ_ij Delta_i, which does not depend on w. So Ehat_ij(ψ, w) = [δ_ij Delta_i - ⟨P^^{(i)}(-ψ), P^^{(j)}(-w)⟩]/(ψ + w), which is the Givental edge in Delta-normalised form. Since η is symmetric, this gives the **edge symmetry** Ehat_ij(ψ, w) = Ehat_ji(w, ψ) for all i, j, so T(Gamma, c) does not depend on the orientation of the edges. (1.1) is a definition, and the conifold Lemma 1.1 (LSP) is the rank-free algebra behind it; see Lemma 1.4 for the code.'''))
PA.append((162, 171, '**Lemma 1.4', r'''**Lemma 1.4 (code bookkeeping).** [P, by reading hodge_vertex/hv_lib3.m lines 57-60 (hodgeClass), 72-73 (hodgeNormU, usubOf), 75-78 (hCoeffs), 119-130 (vertexOptionsH), 131-138 (prepareGraphH), 142-143 (rowsTimes, evalPointRows); hodge_vertex/hv_lib4.m lines 49, 114 (hodgeNormU4, vertexOptionsH4, at genus 4); lp2_lib.m lines 37-64 (evalPoint, graphValue); lp2_explore/lp2_percolour.wls lines 41, 47, 50, 81-84, 105-106, 116-118.] The value computed at frame-B rows by graphValue, and by the separate series loop of lp2_percolour.wls lines 116-118, which has the same form, equals |Aut|^-1 Sum_c T(Gamma, c). (graphs.m vertexOptions, line 46, has pure ψ weights wk and is not the Hodge vertex.)

*Proof.*
- Vertex integral. hodgeClass[g] expands prod_a Sum_j (-1)^j λ_j u_a^j, and with u_a = 1/w_a, w = {E, -ε1, -ε2} (usubOf; lp2_percolour.wls lines 47, 106) this is Omega_w of s.1.5. hodgeNormU[g, a] is Int Omega_w ψ^a. vertexOptionsH produces the entries {hp, lp, dils, wgt} with wgt = hodgeNormU/(legs! nd!), dils the ordered tuples of degrees >= 2, and hp of total degree <= the budget (mixed degree). So at legs = 0 and nd = m_v, wgt is (1/m_v!) Int Omega_w ψ^a ψ^d as in (0.1).
- The code's vertex factor contains sqD^(2g_v - 2 + n_v + legs + #dil), with legs = 0 (G_g has no external markings; the code's gr["m"] is the number of legs, not the m_v of (0.1)) and sqD = Delta^{1/2}. This is Delta^{g_v-1}, times one sqD per half-edge, times one sqD per dilaton.
- One sqD per dilaton turns Dil[i][n] = (-1)^n Pnum[0,i][[n]] into Dil^_i(n).
- One sqD per half-edge turns Ed into the normalised edge. Indeed pre(m, n) = δ_ij δ_m0 δ_n0 - (-1)^{m+n} P_m η P_n has generating function N_ij(ψ, w) = δ_ij - ⟨P(-ψ), P(-w)⟩, and Ed(m, n) = Sum_{k=0}^{n} (-1)^k pre(m+k+1, n-k) has generating function [N_ij(ψ, w) - N_ij(-w, w)]/(ψ + w). The telescoping is ψ^M - (-w)^M = (ψ + w) Sum_{k<M} (-1)^k ψ^{M-1-k} w^k, and M = 0 drops out, so the Kronecker term never enters. Multiplying by sqD_i sqD_j gives (1.1).
- The frame-B twist is rowsTimes[hCoeffs[wB, -1, KMAX], KMAX] (lp2_percolour.wls line 50), i.e. multiplication by h_w^{-1}.
- The Leaf factors Lf occur only with markings, so not in G_g. ∎

This lemma relates (0.1) to the code. The GW meaning of Dil^ (Givental's T(z) = z(1 - R^{-1}(z) 1)) is **not** derived here; it belongs to Proposition 13.2_P2. Theorems LC_P2 and R_P2 are statements about the decoration sum as defined.'''))
PA.append((173, 173, '**Remark (signs).**', r'''**Remark (signs).** [P] In the normalised form (0.1) no square root occurs. Delta_i = 3 ε1 ε2 E/v_i^2, and P^_a = Delta_i^{1/2} P_a^B does not depend on the sign of Delta_i^{1/2} because P_0^B carries Delta_i^{-1/2}, so that P^_a is built from S, θ, Y_i, C1i, d2 and θ log Delta_i, all single valued. Hence each T(Gamma, c) is well defined on any cover. The compatibility of branches in the code is built in, since sqDsym is transported with the rows. (The sign of c1 is a different matter: see Lemma G.)'''))

# s.1.5
PA.append((180, 180, '(This is the convention', r'''(This is the convention of LSP Prop 13.2, of the code hv_lib3.m lines 2-3 and 57-60, and of the conifold's hvc_lib.m lines 3, 139. The scale law of s.7.4 depends on it: with the unnormalised product the exponent there would change.) It is a pull-back under forgetful maps, so the string and dilaton equations hold for h >= 2.'''))
PA.append((182, 192, 'Definitions, for h >= 2:', r'''Definitions, for h >= 2. Let 𝐑 := C[[t_c]] (or Q(ε1, ε2)[[t_c]], or any t_c-adically complete Q-algebra of coefficients).

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

At colour 0, q_t is in 𝔇 by Lemma D0(iv). Edges do not factorise in general: a loop, or an edge between two rescaled vertices, is a tensor. (1.2) is multilinear in the half-edge slots and is applied term by term, leaf by leaf in the expansions of E00(i) and (M-b), and coefficient by coefficient in the edge variables.'''))
