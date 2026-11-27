# Round-5 patch, part D: s.7.
PD = []
PD.append((590, 592, '**(H0) Dilaton homogeneity.**', r'''**(H0) Dilaton homogeneity on 𝔇.** [P; this makes precise LSP s.13.1(i), which does not state a domain.] Let 𝔇 and 𝔡 be as in s.1.5. For q in 𝔇 and λ in 1 + t_c 𝐑 independent of ψ, one has λq in 𝔇 and

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

Step 3 (integration). Put μ := λ - 1, in t_c 𝐑. Let ϑ be a formal variable, q_ϑ := (1 + ϑμ) q, and g(ϑ) := d^n Fpot_h(q_ϑ)[f]. Each t_c-coefficient of g is polynomial in ϑ. The Taylor definition gives dg/dϑ = μ d^{n+1} Fpot_h(q_ϑ)[f, q], and Step 2 at q_ϑ gives (1 + ϑμ) dg/dϑ = μ (2 - 2h - n) g. The solution with value g(0) at ϑ = 0 is unique, so g(ϑ) = (1 + ϑμ)^{2-2h-n} g(0). Put ϑ = 1. ∎'''))

PD.append((622, 622, 'The stability count', r'''The stability count concerns half-edges only: for G_g there are no external markings, and the "legs" of step (ii) below appear only after the graph has been fixed. The decorations agree with the code: hv_lib3.m vertexOptionsH (lines 119-130) has dilaton degrees >= 2 and dim >= 2 nd, so a genus-0 vertex with m = 0 needs at least 3 half-edges. (Round-1 note N6: the single-vertex graph does need h = g.)'''))

PD.append((634, 643, '(b) *Vanishing.*', r'''(b) *Vanishing.* [C: hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md, Theorem L and Corollaries L1-L3, s.7.5 there, for all h >= 2; there P modulo (MV) and textbook facts; status below] For every h >= 2 and every non-empty multiset R of integers >= 2 with D(R) <= 2h - 3: Xi_h(R) = 0 in Q(ε1, ε2), and at every specialisation with ε1 ε2 E != 0, the diagonal ε1 = ε2 (τ = -1/2) included. Combined with (S_P2) and (C_P2), the same vanishing holds for every non-empty R with parts >= 1 in the range. R must not contain a part 0: Xi(0) = ((2-2h)/E^2) Z[h,0], which is a non-zero element of Q(ε1, ε2) by Theorem VI.

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
- **R1 is not admissible** under the all-orders ruling Q5: its base h <= 7 and its zero count are computer-assisted. The computations at h = 3, 4 (PSM s.1.9) are likewise unnecessary.'''))

PD.append((651, 651, 'Fix Γ and a colouring c.', r'''Fix Γ and a colouring c. Call a colour-0 vertex with g_v >= 2 **rescaled**. For a graph with an edge every vertex has at least one half-edge. Fix ρ as in s.1.2. Dilaton decorations at every vertex are summed first, as in (1.2). Everything below is multilinear in the half-edge slots: for an edge that does not factorise (a loop, or an edge between two rescaled vertices), the statements are applied term by term in the expansion of E00(i) or E00(iii), and they are read coefficientwise in the edge variables of the far end.'''))
PD.append((653, 653, '**(i) Homogeneity.**', r'''**(i) Homogeneity.** [P] At a rescaled vertex v of genus h with half-edge series f_1, ..., f_n, by (1.2), (7.1) with λ = ρ (q_t lies in 𝔇 by Lemma D0(iv), and ρ q_t lies in 𝔇 although its ψ-coefficient is -ρ), and Lemma 1.2(b), the vertex factor is'''))
PD.append((671, 671, 'Taylor-expand', r'''By the definition of Fpot_h on 𝔇 (s.1.5), Fpot_h(q_0 + δ) with δ = Sum_{k>=1} t_c^k q^(k) is its Taylor series at q_0. The expansion is finite at each order, because δ = O(t_c) and each Xi is a finite sum. So the coefficient of t_c^{2-2h+k} of the rescaled vertex factor is a finite combination of Xi_h(R) with all parts j >= 1 and D(R) <= k. There are no λ_0 terms and no ψ-polynomial remainders (s.7.2). Each half-edge contributes at least one leaf, so r >= n_v >= 1.'''))
PD.append((687, 687, '**Status of the conclusion.**', r'''**Status of the conclusion.** Theorem R_P2 holds for the decoration sum as defined (s.1.4), per colouring, for every genus, modulo Lemma A_P2(b) for 2 <= h <= g-1. Genus 2 uses no Lemma A at all. No BS input enters (s.0.3). Identification with GW needs (C1)-(C5) and Proposition 13.2_P2, which remains open.'''))

PD.append((691, 703, '**Vertex identity (VI).**', r'''**Theorem VI (vertex identity (VI)).** [P mod LLZ eqs. (1), (5), the boson-fermion correspondence, string/dilaton; proof in s.7.8] For g >= 2,

    (-ε1 ε2 E^3)^{g-1} Z[g, 0] = β_{2g-2}(ε1, ε2),     β_n := (n-1)! [t^n] 1/(4 sinh(ε1 t/2) sinh(ε2 t/2)).

At ε2 = 1, ε1 = x, with c = 1 + x = E, this reads (-x c^3)^{g-1} Z[g, 0] = gap_g, where gap_g = (2g-3)! [t^{2g}] t^2/(4 sinh(xt/2) sinh(t/2)) = (2g-3)! P_g(x,1)/x. This is conifold F4 at n = 0.

Z[g, 0] is a pure Hodge integral: it depends only on Omega_w and the scale E. The proof in s.7.8 uses no TR, no Eynard/ABDKS, no Lemma A, no Teleman and no finite-order input.

Cross-checks, not proof steps:
- (B-a) Lemma A round 2, route R2 Corollary 7 (gap formula, lemmaA/round2/round2_status.md s.1.2), by Hodge-side TR. It is equivalent to Theorem VI by the short conversion of s.7.8, Part A. Its ABDKS/Eynard normalisation was not derived for all g, so it is demoted to a cross-check.
- (B-c) conifold Proposition 13.2 + BS prop:ResConifold (main.tex line 1708; "Prop. 3.1" of arXiv:2410.00118) (H) + Corollary R with Lemma A(<= g-1). This is TR-free, but the Teleman input of conifold Prop 13.2 was read through summaries only (conifold critic gap list), and it keeps a BS citation. So it is a consistency remark. See also BS prop:cfgraphsum (statement at 5114, eq:betagiv 5122, label 6366), which is related and which has not been reconciled with the single-vertex value Z[g,0].

The recorded proof is s.7.8.'''))

PD.append((705, 705, '**Proof of Corollary SV_P2.**', r'''**Proof of Corollary SV_P2.** [P mod M_irr(h), h <= g, and Theorem VI]'''))

PD.append((723, 723, '- **(b) Colouring sums.**', r'''- **(b) Colouring sums.** None are used; every step is per colouring. Each T(Γ, c) contains no square root (s.1.4).'''))
PD.append((724, 724, '- **(c) Per-colouring', r'''- **(c) Per-colouring statement.** Only the label of the singular colour depends on the branch. Near y = -1/27 the roots t_i are analytic and distinct (Lemma 2.1), so labels are locally constant. Continuation is along the real branch y in (-1/27, 0) (round-1 N8, journal wf_e928f771-9e3).'''))
PD.append((727, 727, '- **(f) ε1, ε2', r'''- **(f) ε1, ε2 indeterminate or complex.** Theorem A and the proof of Proposition 2.10 use only E != 0 (and a square root of 3 ε1 ε2 E in S0), so they work verbatim over C(t) for every complex ε with ε1 ε2 E != 0. The statements of s.3-s.6 follow algebraically. Lemma A and Theorem VI hold identically in τ and specialise, through the Laurent polynomial of s.7.4(b) and through (7.2). But the non-vanishing of Z[h,0] does not specialise: at h = 2, Z_2(τ) = -(24 + 28s + 7s^2)/5760 with s = 1/σ_τ, by Theorem VI, and this vanishes at two real values σ_τ < -1/4, i.e. at complex τ. No proof step uses the non-vanishing.'''))

PD.append((730, 730, '- **(i) Theorem R as stated', r'''- **(i) Theorem R as stated in LSP/PSM** assumes M_irr(h) up to h = g; for edge graphs only h <= g - 1 is needed (Lemma 7.1). The documentation of the conifold should be corrected accordingly (round2_status.md A4).

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

Status. One referee pass (logic, literature and independent computation in one lens: sound with minor fixes, all applied above). The referee confirmed the value by an MV-free computation (Mumford/Faber recursion and DVV) at h <= 4 at three τ and at h = 5 at one τ [E]. Scope: Step 1 settles (B1) for all g and μ, the same derivation as Proposition MV(A) of the MVmatch report; it does not by itself establish Lemma A_P2(b), whose status is given in s.7.4.'''))
