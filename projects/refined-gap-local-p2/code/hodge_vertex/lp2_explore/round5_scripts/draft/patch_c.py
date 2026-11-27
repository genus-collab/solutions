# Round-5 patch, part C: s.3-s.6.
PC = []
PC.append((373, 373, 'Put σ := ρ', r'''Put σ := ρ P^^{(0)}(z = -ψ) and σ^(k) := [t_c^k] σ. By Corollary 2.8 (the colour-0, frame-B case), the σ^(k) lie in the space V := C[ψ] ⊕ span{λ_j : j >= 1}; partial fractions make the decomposition unique. Write σ_0^(k) = Sum_j β_j^(k) λ_j. Translating Corollary 2.9 and Proposition 2.10:'''))
PC.append((379, 379, '- **(K=1)**', r'''- **(K=1)** With ' meaning ψ -> -ψ: η^00 σ_0 σ_0' + η^11 σ_1 σ_1' + η^12 (σ_1 σ_2' + σ_2 σ_1') = 3 ε1 ε2 E^3/t_c^2 as an identity of Laurent expansions in t_c with coefficients in V ⊗ V. This is Proposition 2.10 at i = j = 0 and Lemma 1.2(b). Since Proposition 2.10 is an identity of rational functions of t at each z-order, no continuation or branch argument is needed; Corollary 2.8 gives the coefficients in V.'''))
PC.append((383, 383, '**Lemma D0.**', r'''**Lemma D0.** [P, from Corollary 2.9]'''))
PC.append((398, 398, '### 3.3 Closure', r'''### 3.3 Closure (conifold-frame QDE). For general c1, not used in the proof of Theorems LC_P2, R_P2.

The c1 = 1 form of the closure is the first-order system S1 of Proposition 2.10, which is derived there directly from Theorem A. So Proposition 3.1 itself is not an input.'''))
PC.append((400, 408, '**Proposition 3.1', r'''**Proposition 3.1 (round 1 closure; independently re-derived by the round-1 verifier and, for general c1 including the z^1 cancellation and self-adjointness, by the round-5 modrel author).** [P, from Theorem A] Let c1 be any invertible section, XX = θ log c1, and 𝐫 = (P_0^B, P_1^B, P_2^B) the vector of unnormalised rows of Definition 1.3 at any colour (the z-only factor h_w does not matter). With Th := Y + zθ,

    C1i Th 𝐫 = (A_0(y, C1i) - z^2 c1 Rd E_21) 𝐫,

where E_21 is the elementary matrix with entry (row 2, column 1), counting rows and columns from 0, and where
- A_0 = [[0, 1, 0], [0, C1i d2, C1i^3/Δ_c], [-E^3/27, -c1(E d2 + Δ_c d2^2 + E^2/3), -C1i(E/Δ_c + d2)]] is z-free, η-self-adjoint and independent of XX;
- Rd := Δ_c θXX + Δ_c XX^2 + 27y XX + 6y is the Riccati defect. c1 solves L2b iff Rd = 0.

So at c1 = C1cf the conifold-frame QDE z d/dt_c 𝐫 = (A_0 - C1i Y_i) 𝐫 holds at colour i, with A_0 z-free and self-adjoint.'''))
PC.append((412, 412, '**On-shell.**', r'''**On-shell.** XX on the conifold section is the unique Riccati solution analytic at Δ_c = 0, XX(0) = -2/9, because the recursion (n+1) x_n = n x_{n-1} + [XX^2]_{n-1} + (2/9) δ_{n1} is non-resonant. Off shell, the z-even self-adjoint term E_21 drops out of the pairing (which is why symplecticity holds off shell, Proposition 2.10), but it produces the pole 3 ε1 ε2 E^2 (XX(0) + 2/9)/t_c of Ehat_00 [E]. (Round 1's ε1 ε2 E^2 (XX(0) + 2/9) is per 1/s, with t_c ~ 3s; it is the same statement.) The proof below does not use closure for general c1, or the flow. On-shell enters only through c1 = θ t_c, i.e. (I1), (I2) of Lemma 4.2, combined with Z-infinity; it is measured by H2(d).'''))
PC.append((415, 415, '## 4. Lemma H2', r'''## 4. Proposition H2'''))
PC.append((428, 429, '- On shell', r'''- On shell (c1 = g), with D_c := Ŷ + z/t_c + z d/dt_c:
  - (I1) σ_1 = D_c σ_0, where Ŷ := C1i Y = E d log v/dt_c = E/t_c + E/18 + O(t_c);'''))
PC.append((432, 432, 'Remark (not an input)', r'''Remark (not an input). This agrees with round 1's L_1 (here Ŷ, renamed to avoid a clash with the leaf statement L1), κ, μ.'''))
PC.append((436, 438, '    [t_c^n] D_c F', r'''    [t_c^n] D_c F = (E - (n+2)ψ) F^(n+1) + Sum_{m>=0} Ŷ^(m) F^(n-m),

and (LAM2) holds. Here Ŷ^(m) := [t_c^m] (Ŷ - E/t_c).'''))
PC.append((440, 440, '**Proposition H2.**', r'''**Proposition H2.** [P; the inputs of s.3.1 are all P]'''))
PC.append((451, 451, '- n >= 0: by (LAM2)', r'''- n >= 0: by (LAM2) (with N = n+2), (E - (n+2)ψ) σ_0^(n+1) is the constant (n+2) Sum_j β_j^(n+1)/j, which vanishes by Z-infinity since n+1 >= 1, plus terms in Λ_{n+1} (the top leaf λ_{n+2} cancels because (E - (n+2)ψ) λ_{n+2} = 1). The Ŷ terms lie in Λ_{n-m+1}, and the m = n+1 term vanishes because σ_0^(-1) = 0. The E/t_c term of Ŷ is built into the prefactor E - (n+2)ψ.'''))
PC.append((490, 490, '**Theorem E00', r'''**Theorem E00 (R3:Theorem E00).** [P, using s.3.1 and H2 Steps 1-2] By (1.1) and K = 1, Ehat_00 = [Delta_0 - ⟨τ(ψ), τ(ψ')⟩]/(ψ + ψ'), τ(ψ) := P^(-ψ). Put 𝒴 := t_c σ and'''))
PC.append((504, 504, '- The hypotheses (H1)-(H3)', r'''- The hypotheses (H1)-(H3) of round 1 map as follows: H1 = Definition 1.3 (the identification with BS's rows, Theorem 2.7, is no longer needed here); H2 = Proposition H2; H3 = K = 1 (Proposition 2.10).'''))
PC.append((512, 512, 'Throughout, j, j', r'''Throughout, j, j' in {1, 2}, and v_{j,k} := [t_c^k] P^^{(j)}, a vector of formal power series in w (at z = -w). It is not a vector of rational functions: since v_j(0) != 0, the l-sum in S(v_j, -w) does not terminate (Corollary 2.8). Since Delta_j^{1/2} is analytic and non-zero at t_c = 0, Corollary 2.9(R) gives v_{j,k} = 0 for k < 0. (M-a)-(M-c) are taken [w^b]-wise, so they are not affected. Moreover, the component 2 of v_{j,0} is 0 identically: in Definition 1.3, P_2 has the factor Δ_c and the bracket is analytic at c1 = C1cf.'''))
PC.append((514, 514, '**Theorem M', r'''**Theorem M (R3:Theorem M).** [P, using Corollary 2.9(R), Proposition 2.10, H2]'''))
PC.append((520, 523, '  At t_c = 0 this is rank one', r'''  At t_c = 0 this is rank one, λ_1(ψ) χ_j(w), with χ_j(w) = -⟨(E, E^2/9, 0), v_{j,0}(-w)⟩/(E + w) (from σ^(0), H2(c)).
- (M-c) For every integer e: [t_c^k][w^b](ρ^e Ehat_0j) lies in Λ_{k+1}(ψ), with no constant.
- (M-d) Ehat_j0(ψ, w) = Ehat_0j(w, ψ).
- (M-e) Ehat_jj', Dil^_j and Delta_j are analytic at the conifold; Delta_j(0) = ε1 ε2 E^3 ω^{2j}/(3(1-ω^j)^2) != 0; and Ehat_jj'(ψ, w) = Ehat_j'j(w, ψ).'''))
PC.append((531, 533, '(M-d). Proposition 2.10', r'''(M-d). This is the edge symmetry of s.1.4, from Proposition 2.10 at (i, j) = (0, j). That is ⟨P^^{(0)}(u), P^^{(j)}(-u)⟩ = 0, as an identity of Laurent expansions in t_c: by the drop-out, the pairing depends only on P_0 and each z-coefficient is rational in t. So P_0j(-w, w) = 0 = P_j0(-w, w), and η is symmetric.

(M-e). Each Ehat_jj'(m, n) is a polynomial in finitely many row coefficients, analytic by (R). Dil^_j uses component 0 only. The value of Delta_j uses 1 + 3t_j -> 1 - ω^j != 0, by s.1.1. The symmetry is the edge symmetry of s.1.4. ∎'''))
PC.append((538, 538, '- The rank-one form', r'''- The rank-one form at t_c = 0 was also verified directly in all entries [E, round-3 referee script R4], and the formula for χ_j reproduces the stored χ_1, χ_2 at (7,3) [E, round-4 critic 2]. Round 3's "conjectural rank-2 relation" and the inconsistent third entry of its stored vector are explained by a one-liner: the third entry multiplies a zero row. Indeed, v_{j,0} has component 2 equal to 0 by Definition 1.3, and the order-t_c^{-1} part of orthogonality with σ^(-1) = (0, E, -E^2) gives the same: ⟨σ^(-1), v_{j,0}⟩ = -3 ε1 ε2 E (v_{j,0})_2.'''))
PC.append((540, 540, '**Corollary 6.1', r'''**Corollary 6.1 (edge inputs).** [P] Theorem LC_P2 holds, with no BS input (Proposition 2.10 is proved from Theorem A).'''))
PC.append((557, 557, '    Ehat_0j    λ_1(ψ) χ_j(w)', r'''    Ehat_0j    λ_1(ψ) χ_j(w), rank one, χ_j(w) = -⟨(E, E^2/9, 0), v_{j,0}(-w)⟩/(E + w)'''))
