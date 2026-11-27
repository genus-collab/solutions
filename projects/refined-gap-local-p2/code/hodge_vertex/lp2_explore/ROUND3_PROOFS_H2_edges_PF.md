

# Round 3 consolidated write-up: PF and F1-F5, Lemma H2, the residue lemma, mixed edges

Labels: [P] proved here; [R2:x] round-2 result x of ROUND2_PROOFS_rows_on_TF.md; [BS:n] BS v2 main.tex at about line n. Exact computer checks are evidence only and are never used as proof steps.

## 0. Setting and conventions

- **The torus and the colours.** We work on T_F with E = e1 + e2, e1·e2·E ≠ 0 and α = E/3. Colour i has t_i = ω^i y^{1/3}. Put v_i = (1+3t_i)/(E t_i), V = E v_i, Y_i = 1/v_i - E/3 and θ u_i = Y_i, where θ = y d/dy.
- **Normalised rows.** P̂_i := Δ_i^{1/2} P_i^B, with Δ_i^{1/2} = √(3e1e2E)/v_i [R2:s.0, Prop 8 step 3]. In particular Δ_0 = 3e1e2E³t²/(1+3t)².
- **The pairing.** ⟨a,b⟩ = Σ η^{ab} a_a b_b, with η^00 = e1e2E³/9, η^11 = -3Ee1e2, η^12 = η^21 = -3e1e2, and all other entries 0 (lp2_lib etainvF). Note η^11 = E·η^12.
- **Leaves.** λ_j(ψ) = 1/(E - jψ) for j ∈ ℤ∖{0}. Λ_m = span{λ_1, …, λ_m}, Λ_0 = 0 and Λ_m^+ = ℂ + Λ_m. For the bounds below we read Λ_m as Λ_{max(m,0)}. The edge variable is ψ = -z.
- **Conifold coordinate.** Δ_c = 1 + 27y. We normalise t_c = Δ_c + O(Δ_c²), so that C1cf = θ t_c and C1cf(0) = -1 [R2:s.0; lp2_lib C1cf = 27y dt_c/dΔ_c]. On the real branch v_0 = t_c·u(t_c), with u a z-free unit and u(0) = -1/E.
- **The unit ρ.** ρ := -Ev_0/t_c. It is analytic in t_c with ρ(0) = 1 (this uses the normalisation of t_c), and ρ²Δ_0 = 3e1e2E³/t_c². Multiplying by ρ^e (e ∈ ℤ) preserves analyticity and the filtration by Λ_{k+1}.
- **Polarisation.** The polarisation is c1 = C1cf with XX = θ log C1cf, the same c1 at every colour; this is legitimate because c1 depends on y only. In the row formulas P1 = C1i·Q1, we have C1i = 1/c1.

**Lemma CL (carrying identities to the conifold).** [P]
- Each z^n coefficient of any row at any colour is a polynomial in v_i, multiplied by functions of t_c that are analytic at the conifold and z-free (C1i, θv, Δ_c/t_c, d2, Δ_i^{±1/2}). So it is a meromorphic germ at y = -1/27. An identity between such series that holds on a punctured neighbourhood therefore holds as an identity of Laurent expansions in t_c.
- Since v_0 = O(t_c), each [t_c^k] coefficient involves only finitely many v-orders. So [z^n] and [t_c^k] commute, and each [t_c^k] coefficient is a rational function of z.

## Part I. Picard-Fuchs for the closed form, and formal stationary phase

**Definitions.**
- (Rec): r_{-1} = 0, r_0 = 1, and 27(m+1) r_{m+1} = 3(3m²+3m+1) r_m - m³ r_{m-1} for m ≥ 0.
- S(v,z) := Σ_l r_l (vz)^l / ∏_{j=1}^{l+1} (1 + jz/E).
- With a := 1 + E/z and (a)_n the Pochhammer symbol, put e_n := V^n/(a)_n for n ∈ ℤ. Then vS = z^{-1} Σ_{n≥1} r_{n-1} e_n.

**Theorem A.** [P] Let P := e^{-u_i/z} L_pf e^{u_i/z}, with L_pf as in [BS:2627-2632]. Let σ e_n = e_{n-1}, N e_n = n·e_n and τ = 3σ - N. Then:
- P = (z³/27)(τ³ - 27σ³) as operators.
- P(Σ b_n e_n) = -(z³/27) Σ_m [27(m+1) b_{m+2} - 3(3m²+3m+1) b_{m+1} + m³ b_m] e_m.
- Hence Δ_i^{-1/2} h_w(z) S(v_i,z) e^{u_i/z} solves eq:pf identically, at i = 0, 1, 2.

*Proof.*
- (A0) Conjugation gives P = T³ + 27y ∏_{j=0}^{2} (T - E/3 + jz/3), where T = zθ + Y + α.
- (A1) Three facts hold at every colour: θV = 1 - V/3, Y + α = E/V, and y(V-3)³ = 1.
- (A2) We have V e_n = (a+n) e_{n+1} and z(a+n-1) = E + nz. Hence T e_n = z e_{n-1} - (nz/3) e_n, that is, T = (z/3)τ, and T - E/3 + jz/3 = -(Ω - jz)/3 with Ω := E + zN - 3zσ.
- (A3) Let M be multiplication by z(V-3). Then Ω = Mσ and [σ, M] = z. It follows that M²σ² = Ω(Ω - z) and M³σ³ = Ω(Ω - z)(Ω - 2z). Since yM³ = z³, we get P = (z³/27)(τ³ - 27σ³).
- (A4) Expanding, τ³ e_n = 27e_{n-3} - 27(n-1)e_{n-2} + 3(3n²-3n+1)e_{n-1} - n³e_n, which gives the displayed coefficient. With b_n = r_{n-1}, that coefficient vanishes by (Rec) for m ≥ 0, and trivially for m ≤ -1.
- (A4', z-adic rigour) Since 1/(a)_n = ∏_{j≤n} z/(E + jz), each e_n with n ≥ 0 lies in V^n z^n Q[[z]], so Σ b_n e_n converges z-adically. P has coefficients polynomial in z and never lowers the z-order. So [z^k] P(vS) is a finite expression, a Laurent polynomial in V and (V-3), and it vanishes in Q(t_i).
- Finally, Δ_i^{-1/2}, h_w(z) and the factor z^{-1} commute with P. ∎

**Corollary A5 (uniqueness).** [P] Suppose b_n ∈ ℂ((z)), b_n = 0 for n ≤ 0, and P(Σ b_n e_n) = 0. Then b_n = b_1(z) r_{n-1}, so Σ b_n e_n = z·b_1(z)·vS. PF alone therefore fixes the row only up to a factor depending on z alone; c(z) = 1 comes from [R2:Prop 8] together with the transport from y → 0.

**Part C (formal stationary phase).** [P]
- *The ring.* Let A be the ring of sums Σ c_{k,β} z^k ξ^β with finitely many terms having k + |β|/2 ≤ N, for each N. Put ⟨f⟩ := e^Δ f |_{ξ=0}, with Δ = (z/2) H^{ab} ∂_a ∂_b. This is continuous on A. Define Λ(γ) := (det H)^{-1/2} ⟨e^{-φ^ng/z} γ⟩.
- *G1 (Gaussian integration by parts).* ⟨(Hξ)_a f⟩ = z⟨∂_a f⟩, because [Δ, ξ^b] = z H^{bc} ∂_c commutes with Δ.
- *F2 (total derivatives vanish).* Λ(∂_a B - z^{-1}(∂_a φ)B) = 0, by G1.
- *FP (fixed critical point).* Suppose the critical point is at ξ = 0 for every p, with critical value independent of p. Then d/dp Λ_p(γ_p) = Λ_p(γ_p' - z^{-1} φ_p' γ_p). Proof: apply G1 twice to reproduce -(1/2) tr(H^{-1}H') and Δ'. This is termwise legitimate by [R2:Lemma 2].
- *F3 (moving critical point).* Translate the critical point to 0 and apply FP; the x0' terms die by F2. Result: d_p[e^{-φ_0/z} Λ] = e^{-φ_0/z} Λ(∂_p γ - z^{-1}(∂_p φ)γ).
- *F1 (coordinate invariance).* With the branch convention (det H')^{1/2} = det Dg(0)·(det H)^{1/2}:
  - The linear case is direct.
  - In general, interpolate by g_s(ξ) = x0 + s^{-1}(g(sξ) - x0). Then H_s is constant in s, and X_s = (Dg_s)^{-1} ∂_s g_s vanishes at 0. One has φ_s' = X(φ_s) and γ_s' = X(γ_s) + div(X)γ_s.
  - By FP and F2, Φ'(s) = 0. Each z-coefficient is polynomial in s, so Φ(1) = Φ(0).
- *F4 (Fubini).* If the phase splits exactly, Δ splits into commuting parts. The amplitude need not split.
- *F5 (formal Stirling).* For α ≠ 0, Asym_{q=α}[e^{-(q - α log q)/z} q^{β-1} dq] = α^{β-1/2} Θ(z/α, β), with the branch fixed by the F1 convention.
  - Scaling gives the form α^{β-1/2} F(w,β) with w = z/α.
  - F2 with B = q^β gives the shift relation (S): F(β+1) = (1 + βw)F(β).
  - F3 in α gives the derivative relation (D): ∂_β F = (β - 1/2)wF - w²∂_w F.
  - Write log F = Σ g_k w^k. (S) makes the difference of two solutions constant in β, and (D) at index k+1 kills those constants, so the solution is unique.
  - Θ satisfies (S) and (D), using B'_{k+1} = (k+1)B_k. Therefore F = Θ. This is the formal version of [BS:exmpl:stirling].

**Corollary B.** BS's Wick numbers [BS:5252-5262] equal the r_l of (Rec) for all l. *Proof:*
- (i) e^{u/z}A_0 solves eq:pf [R2:Lemma 3]. Its step 2 is F3 applied three times, with amplitudes that are θ-derivatives, plus F2 on θ-images of total derivatives; θ commutes with ∂_x. [P mod writing out]
- (ii) A_0 = K_0·Θ(-z/E)·vS_Wick [R2:Prop 8 steps 1-2]. This uses the s-chart W = s1³ - 3s1s2 + (3 + y^{-1/3})s3 - E log s3, F1, F4, and F5 with α = -E, β = -l-1, where Θ(w,-m) = Θ(w)/∏_{j≤m}(1 - jw).
  - The (-1)^{l+1} from s3 = -q/V cancels against (-E)^{-l-1}, so the prefactor has no l-dependent sign.
  - The branch signs are absorbed into K_0.
- (iii) Θ(-z/E) is invertible, so P(vS_Wick) = 0. Corollary A5 with r^W_0 = 1 gives r^W_l = r_l. ∎

**Consequence.** [R2:Lemma 3, Prop 5, Theorem 6, Prop 8, Cor 9, Prop 10] now hold modulo only BS C1-C3 and step (i) above.

## Part II. Lemma H2

Notation at colour 0: v = v_0. Then θv = 1/E - v/3, Y = -E/3 + 1/v, θ log v = Y/E, and Δ_0 ∝ v^{-2}. Put P̂_0 = S [R2:Prop 8], σ := ρ·P̂(z = -ψ), and σ^(k) := [t_c^k]σ. The coefficients σ^(k) lie in V := ℂ[ψ] ⊕ span{λ_j : j ≥ 1}, where partial fractions make the decomposition unique (Lemma CL).

**Lemma 1 (closed forms for rows 1 and 2, any section c1).** [P] Put S^[m] := Σ_l r_l (vz)^l / ∏_{j=1}^{l+1-m} (1 + jz/E). Then
- P̂_1 = E·C1i·(θv/v)·S^[1];
- P̂_2 = E·c1·Δ_c·θv·[E θv S^[2]/v² + (z θ log(C1i θv) - d2) S^[1]/v].

*Proof.* Use the identity Y + zθ log v = (E+z)θv/v, together with the termwise identities (E + (l+1)z)/∏_{j≤l+1} = E/∏_{j≤l} and (E + lz)/∏_{j≤l} = E/∏_{j≤l-1}. This applies to the rows of [R2:Definition 7]. ∎

**Lemma 2.** [P] Put g = θt_c and D_g := Y + z·g·(1/t_c + d/dt_c). Then:
- (G) For any section: σ_1 = C1i·D_g σ_0 and σ_2 = Δ_c c1² (D_g - d2) σ_1. This is conjugation by ρ, using θ log ρ = θ log v - θ log t_c and θ = g d/dt_c.
- On shell (c1 = g), with D_c := L_1 + z/t_c + z d/dt_c:
  - (I1) σ_1 = D_c σ_0, where L_1 = C1i·Y = E d log v/dt_c = E/t_c + E/18 + O(t_c);
  - (I2) σ_2 = κ D_c σ_1 + μ σ_1, where κ = c1³Δ_c = -t_c + O(t_c²) and μ = (E/2)(c1² - 1) = O(t_c).
- *Remark (not an input):* this agrees with round 1's L_1, κ and μ.

**Lemma 3.** [P] For F ∈ t_c^{-1} V[[t_c]]:
- [t_c^n] D_c F = (E - (n+2)ψ) F^(n+1) + Σ_{m≥0} L_1^(m) F^(n-m).
- (LAM): (E - Nψ)λ_j = N/j + (1 - N/j) E λ_j. In particular (E - Nψ)λ_N = 1.

**Inputs.** All are [R2], made unconditional by Part I and extended to the conifold by Lemma CL:
- (L1) σ_0^(-1) = 0 and σ_0^(k) ∈ Λ_{k+1};
- (F2) σ_0^(0) = Eλ_1;
- (Z∞) writing σ_0^(k) = Σ_j β_j^(k) λ_j, one has Σ_j β_j^(k)/j = E δ_{k0};
- (K=1) η^00 s0 s0' + η^11 s1 s1' + η^12 (s1 s2' + s2 s1') = 3e1e2E³/t_c², where ' means ψ → -ψ.

**Proposition H2.** [P mod the inputs]
- (a) σ^(k) = 0 for k ≤ -2, and σ^(-1) = (0, E, -E²).
- (b) For k ≥ 0, σ_0^(k) and σ_1^(k) lie in Λ_{k+1}, and σ_2^(k) lies in Λ_k. There is no ψ-polynomial part and no constant.
- (c) σ^(0) = (Eλ_1, (E²/9)λ_1, 0).
- (d) Within the family c1 = C1cf·f(Δ_c), the ψ-linear coefficient of σ_2 is E Δ_c C1cf² θ(f)/t_c. It vanishes iff f is constant. Analyticity of d2 then forces f = ±1.

*Proof.*
- **Step 1 (σ_1).** Apply (I1) and Lemma 3 with F = σ_0.
  - For n = -1: (E - ψ)Eλ_1 = E.
  - For n ≥ 0: by (LAM), (E - (n+2)ψ)σ_0^(n+1) is the constant (n+2)Σ_j β_j^(n+1)/j, which vanishes by Z∞, plus terms in Λ_{n+1}; the top leaf cancels. The L_1 terms lie in Λ_{n-m+1}, and the m = n+1 term vanishes because σ_0^(-1) = 0.
- **Step 2 (σ_2).** Apply Lemma 3 with F = σ_1.
  - For n ≤ -3 the coefficient [t_c^n] D_c σ_1 is 0. It is E² at n = -2 and a constant at n = -1.
  - For n ≥ 0 it lies in Λ^+_{n+1}. A ψ-linear term could only come from (E - (n+2)ψ)·const(σ_1^(n+1)), and that is zero: by Step 1 for n+1 ≥ 0, and because n+2 = 0 for n = -2.
  - By (I2), using κ^(0) = μ^(0) = 0 and κ^(1) = -1: σ_2^(k) ∈ Λ^+_k for k ≥ 0, and σ_2^(-1) = -E².
- **Step 3 (constants).** Each component of σ is bounded as ψ → ∞, with limit equal to its constant part.
  - In (K=1), whose arguments are ψ and -ψ, take the single limit ψ → ∞. Using η^11 = Eη^12 and η^22 = 0: -3e1e2[E(E/t_c)² + 2(E/t_c)σ_2^∞] = 3e1e2E³/t_c².
  - Hence σ_2^∞ = -E²/t_c, and every constant of σ_2^(k) with k ≥ 0 is 0.
- **Step 4 ((c)).** With r_1 = 1/9 and v = -t_c/E - t_c²/(18E) + O(t_c³): σ_0^(1) = -(E/18)λ_1 + (E/9)λ_2. Then σ_1^(0) = (E - 2ψ)σ_0^(1) + (E²/18)λ_1 = (E²/9)λ_1.
- **Step 5 ((d)).** Coefficientwise in v: lim_{z→∞} S^[1] = S_∞'/E and lim_{z→∞} S^[2]/z = v S_∞''/E², with S_∞ = Σ r_l (Ev)^{l+1}/(l+1)!.
  - So lim_{z→∞} P̂_2/z = Δ_c/(C1i² v)·θ(C1i θ S_∞).
  - By Z∞, S_∞ = -t_c, so C1i θ S_∞ = -C1cf/c1. ∎

## Part III. Residue lemma and the singular loop

**Lemma R.** [P] Let f_0 = 1 and f_j = λ_j. Let P(ψ,ψ') = Σ_{i,l=0}^{K} C_il f_i(ψ) f_l(ψ') with C_il = C_li, and suppose κ + P(ψ,-ψ) = 0 for a constant κ. Then:
- κ = -C_00;
- E C_i0 + Σ_{l≥1} i C_il/(i+l) = 0 for every i ≥ 1;
- [κ + P(ψ,ψ')]/(ψ+ψ') = Σ_{i,l≥1} il/(E(i+l)) C_il λ_i(ψ) λ_l(ψ').

*Proof.* Use λ_i(ψ) - λ_i(-ψ') = i(ψ+ψ')λ_i(ψ)λ_{-i}(ψ') and λ_{-i}λ_l = (iλ_{-i} + lλ_l)/(E(i+l)). On the diagonal, P(ψ,-ψ) decomposes in the basis {1, λ_j : j ≠ 0}, which is linearly independent because the poles are distinct. ∎

**Theorem E00.** [P mod H2 Steps 1-2 and K=1] Define
Ê_00 = [Δ_0 - ⟨τ(ψ), τ(ψ')⟩]/(ψ+ψ'), with τ(ψ) = P̂(-ψ).
This is the Givental edge in Δ-normalised form (edge variable ψ = -z; matches edges_sing.pl lines 8-22). Put Y := t_c σ and Ē := t_c²ρ²Ê_00 = [3e1e2E³ - ⟨Y,Y'⟩]/(ψ+ψ'). Then:
- (i) Ē = Σ_{n≥2} t_c^n Σ_{i,l≥1, i+l≤n} il/(E(i+l)) C^(n)_il λ_iλ_l', where C^(n)_il = -Σ_{k1+k2=n-2} Σ_{a,b} η^{ab} [λ_i]σ_a^(k1) [λ_l]σ_b^(k2).
- (ii) Ê_00 lies in Q(e)[[ψ,ψ']][[t_c]]; in particular it is regular at the conifold.
- (iii) [t_c^k](ρ^e Ê_00) lies in span{λ_iλ_l' : i,l ≥ 1, (i-1)+(l-1) ≤ k}, and it is symmetric.
- (iv) Ê_00|_{t_c=0} = -(e1e2E⁴/27) λ_1 λ_1'.
- (v) [t_c^k][ψ'^b](ρ^e Ê_00) lies in Λ_{k+1}(ψ).

*Proof.*
- Apply Lemma R order by order, with κ = 3e1e2E³δ_{n0}. The hypothesis Ň(ψ,-ψ) = 0 is K=1.
- Y^(0) = (0, E, -E²) is constant, and the constants of Y^(a) feed only C_i0, C_0l and C_00. By H2, i ≤ a and l ≤ b with a + b = n, which gives (i).
- At n = 0, 1 there are no admissible pairs, which gives (ii). (iii) and (v) follow by multiplying by the z-free unit ρ^{e-2}.
- (iv): C_11 = -⟨b,b⟩ with b = (E, E²/9, 0) and ⟨b,b⟩ = 2e1e2E⁵/27. ∎

## Part IV. Mixed and regular edges

Put P_ij(z,w) = -⟨P̂_i(-z), P̂_j(-w)⟩ and E_ij = [P_ij(z,w) - P_ij(-w,w)]/(z+w), the conifold Lemma 1.1, which is pure algebra. Write v_{j,k} = [t_c^k]P̂_j; by [R2:Cor 9(R)] it vanishes for k < 0.

**Lemma G (gauge).** [P] With the same c1 at both colours, η^11 P1P1' + η^12(P1P2' + P2P1') depends only on (P0, XX), because η^11 = Eη^12. So every E_ij is the same in the hybrid gauge (c1 = 1, XX = θ log C1cf) and in the true gauge. The map between the two gauges is z-free, analytic and invertible (C1cf(0) = -1), so H2 holds in either gauge.

**Theorem M.** [P mod R2:Cor 9(R), Prop 10, and H2] For j = 1, 2:
- (M-a) ρE_0j and E_0j lie in ℂ[[t_c]][[ψ,w]].
- (M-b) Write σ^(k1) = c^(k1) + Σ_{l=1}^{k1+1} S_l^(k1) λ_l. Then
  [t_c^k](ρE_0j) = -Σ_{k1=0}^{k} Σ_{l=1}^{k1+1} l ⟨S_l^(k1), v_{j,k-k1}(-w)⟩ λ_l(ψ)/(E + lw).
  At t_c = 0 this is rank one, λ_1(ψ)χ_j(w).
- (M-c) For every integer e: [t_c^k][w^b](ρ^e E_0j) lies in Λ_{k+1}, with no constant.
- (M-d) E_j0(z,w) = E_0j(w,z).
- (M-e) E_jj' (j, j' ∈ {1,2}), Dil_j and Δ_j are analytic at the conifold, and Δ_j(0) = e1e2E³ω^{2j}/(3(1-ω^j)²) ≠ 0. Also E_21 = E_12^T.

*Proof.*
- **(M-a), (M-b).** Put N = ρP_0j = -⟨σ(ψ), P̂_j(-w)⟩. Then (ψ+w)ρE_0j = N(ψ,w) - N(-w,w).
  - Every ψ-constant cancels: all of σ^(-1), and the constants allowed in Λ^+.
  - The λ_l part contributes -(ψ+w) l ⟨S_l, v_{j,k2}(-w)⟩ λ_l(ψ)λ_l(-w), using (LAM6).
  - Dividing by ψ+w gives (M-b). The coefficients for k ≤ -1 vanish because σ^(-1) is constant.
  - By Lemma CL, each E_0j(m,n) is a meromorphic germ with nonnegative t_c-exponents, hence analytic.
- **(M-c).** Take [w^b] in (M-b) and use the unit ρ.
- **(M-d).** Prop 10, carried to the conifold by Lemma CL (or directly, since by Lemma G the pairing depends only on P0 and XX), gives ⟨P̂_0(u), P̂_j(-u)⟩ = 0. So P_0j(-z,z) = 0 = P_j0(-w,w), and η is symmetric.
- **(M-e).** Each E_jj'(m,n) is a polynomial in finitely many row coefficients, which are analytic by (R). The value of Δ_j uses 1 + 3t_j → 1 - ω^j ≠ 0.
- *Remark (heuristic, not load-bearing).* A ψ-linear term Aψ in σ would add the ψ-constant -⟨A, P̂_j(-w)⟩ to ρE_0j. This is what the off-shell controls show. ∎

## Part V. Edge inputs for Theorem R_P2

**Corollary.** [P mod BS C1-C3 and Part I (i)] The edge inputs of the port of Theorem R (leaf_structure_proof s.13.1; c → E, xc³ → -3e1e2E³, t → t_c) hold:
- (T-F3) All nine E_ij are regular: E00(ii), (M-a), (M-d), (M-e).
- (T-L2) Legs to the regular colours, in both orientations: (M-c), (M-d).
- (T-L3) Loops, edges between two rescaled vertices, and legs to unrescaled singular vertices: E00(iii), (v).
- Also F1 at all colours: Δ_0 = 3e1e2E³/(t_cρ)², and Δ_j ≠ 0.

The inputs that are not edge statements remain: L1 and F3 for Dil_0 (round-2 Z∞), the Hodge side (dilaton, Theorem C, string identity), Lemma A, and writing out the port of Theorem R itself.

Evidence scripts (not proof steps), under <user home>/AppData/Local/Temp/claude/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/lp2_round3/: `h2/`, `H2_ref2/`, `E0j/`, `E0j_ref2/`, `PF/`, `PF_ref1/`, `PF_ref2/`.