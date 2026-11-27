

# Rows on T_F: definition, continuity, closed form and consequences

Status labels: **[P]** proved here; **[C]** cited (BS v2 label and line); **[S]** standard formal stationary-phase fact; **[E]** finite-order evidence, never a proof step.

## 0. Setting, notation and inputs

**Weights.** The torus weights are ε1, ε2 (called e1, e2 in the code), with E := ε1+ε2 and ε1ε2E ≠ 0. The remaining weights are α = (α0,α1,α2) ∈ C³, and e_k(α) denotes their elementary symmetric functions. The torus T_F is the point α* = (E/3)(1,1,1); on T_F we write α := E/3. Further notation:
- θ := y d/dy, t := y^{1/3}, ω := e^{2πi/3}, Δ_c := 1+27y;
- L̂(Y) := Π_j(Y+α_j) + 27yY³ (eq:LLMM);
- W(Y) := e1(α)Y² + 2e2(α)Y + 3e3(α) (eq:delta1);
- L_pf := Π_j(zθ+α_j) + z³y Π_{j=0}^{2}(3θ+j) (eq:pf, line 2628).

**Colours on T_F.** Put t_i := ω^i t. Then:
- Y_i := −α/(1+3t_i) and q_i := Y_i+α = Et_i/(1+3t_i);
- v_i := 1/q_i = 3/(E+3Y_i);
- W(Y_i) = E q_i², Δ_i := 3ε1ε2W(Y_i) = 3ε1ε2E/v_i², and we set Δ_i^{1/2} := √(3ε1ε2E)/v_i (this is the code's sqDsym).

At colour 0 the conifold t = −1/3 corresponds to v_0 = 0. There v_0 = t_c·(unit) with v_0 ≈ −t_c/E, where t_c = Δ_c + O(Δ_c²), C1cf = θt_c and C1cf(0) = −1.

**Special functions.**
- Θ(x) := exp Σ_{m≥1} B_{2m}x^{2m−1}/(2m(2m−1)), so Θ(−x) = 1/Θ(x).
- h_w(z) := Θ(−z/E)Θ(z/ε1)Θ(z/ε2), which satisfies h_w(z)h_w(−z) = 1.
- N(z) := √(ε1ε2)Θ(−z/ε1)Θ(−z/ε2).
- r_l are BS's Wick numbers (lines 5257–5262): 1, 1/9, 2/81, 14/2187, …
- 𝔖(v,z) := Σ_{l≥0} r_l(vz)^l / Π_{j=1}^{l+1}(1+jz/E).

**Superpotential.** In log coordinates x_j = log q_j (eq:spotP2, line 4301):
- 𝒲 = Σ_j(e^{x_j} − α_jx_j) + y^{−1/3}e^{(x0+x1+x2)/3};
- dΩ = (1/3)d³x;
- P := y^{−1/3}e^{Σx/3}.

Asym denotes BS's normalised formal Laplace expansion (line ~4245).

**Standard facts [S].**
- (F1) Asym of a top form is invariant under local biholomorphisms that fix the critical point.
- (F2) Asym[∂_{x_j}(Be^{𝒲/z})d³x] = 0.
- (F3) θ_y(e^{u/z}Asym[F]) = e^{u/z}Asym[θ_yF], where u = 𝒲(x^cr).

**Cited inputs [C].**
- (C1) For generic α, e^{u_i/z}R^α_{0i} is a formal solution of eq:pf with θu_i = Y_i and Y_i(0) = −α_i. Here R^α_{0i} = S_{0i}e^{−c_i}, and its coefficients continue analytically in the semisimple locus. [eq:Ifun 2621, eq:J=I 2663, eq:asymp 2943 with ~2950, eq:canR 3566]
- (C2) R^α_{0i}|_{y=0} = Δ_i(0)^{−1/2}e^{−c_i}, with e^{−c_i} = Π_{k=0}^{4}Θ(−z/w_i^{(k)}) over the five tangent weights {3α_i, α_j−α_i, α_k−α_i, −ε1, −ε2}. The z^0 coefficient of R^α_{0i} is Ψ_{0i} = Δ_i^{−1/2}. [(Normalisation) 3497, eq:phi 3031, weights 2383]
- (C3) eq:RHb (3573): R_1 = (zθ+Y)R_0/A and R_2 = (zθ+Y−e1(α)E(y))R_1/D. The quantities A, D = X/A², E(y) = (D−X)/2 do not depend on α (eq:ADEFXdef, eq:Ian, lem:modrel 2695). The pairing η is analytic and invertible at α* (eq:etaflat 2419). (Symplecticity) (3496) holds for generic α.
- (C4) prop:algindABC (2854). Used only to make Cf well defined.
- (C5) Givental–Teleman on rational tori (eq:gtH 3799, thm:haegen 3902), lem:lho (3247), and regularity of the y-coefficients of GW_g^{T_max} on the diagonal (line ~4056). Used only in Proposition 11.

## 1. Critical points

**Lemma 1 [P].**
- The critical points of 𝒲 are exactly the points with q_j = α_j+Y and P = −3Y, where L̂(Y) = 0 and Y(α_j+Y) ≠ 0.
- The Hessian satisfies det Hess 𝒲 = W(Y)/3.
- On L̂ = 0 we have W = −Y L̂′(Y).
- At T_F the three roots Y_i are finite, nonzero, simple and nondegenerate if and only if y ∉ {0, −1/27}.
- At y = −1/27 only the colour-0 root escapes, since 1+3ω^it = 1−ω^i ≠ 0 for i = 1, 2.
- θu = Y at each critical point.

*Proof.*
1. ∂_j𝒲 = q_j − α_j + P/3. Put Y := −P/3; then the constraint P³ = Πq/y becomes L̂(Y) = 0.
2. Lifts of Y to log coordinates differ by 2πiZ and have identical jets.
3. Hess = diag(q) − (Y/3)J, with J the all-ones matrix. Using e2(q) = e2 + 2e1Y + 3Y², det = Πq − (Y/3)e2(q) = W/3.
4. From 81yY² = −3Π(Y+α)/Y we get Y L̂′ = Y e2(Y+α) − 3Π(Y+α) = −W.
5. On T_F, (Y+α)/Y = −3t_i, which gives the formulas of §0 and W(Y_i) = Eq_i² ≠ 0.
6. θu = θ𝒲|_{cr} = −P/3 = Y. ∎

## 2. Analyticity of saddle coefficients

**Lemma 2 [P].** Let p range over a parameter domain on which the critical point x_c(p) is nondegenerate and analytic, and choose an analytic branch of (det H)^{−1/2} on a simply connected set. Then every z-coefficient of Asym is analytic in p. For 𝒲, on a simply connected neighbourhood U of the semisimple locus in (y,α), the z^n coefficient ρ_n of 𝒜 := Asym_{x^cr}[e^{𝒲/z}dΩ] lies in W^{−1/2}Q[α,Y][W^{−1}].

*Proof.*
1. The z^n coefficient is the finite sum Σ_{m−k=n, 3k≤2m}(2^m m! k!)^{−1}(H^{ab}∂_a∂_b)^m[(−f_ng)^k g]|_{x_c}, with k ≤ 2n and m ≤ 3n.
2. Each term is a polynomial in the jets at x_c and in adj(H)/det H.
3. For 𝒲 the jets of order |β| = m ≥ 2 are ∂^β𝒲 = q_jδ_{β=me_j} − Y/3^{m−1}, and the amplitude is constant. ∎

## 3. Picard–Fuchs

**Lemma 3 [P mod F2, F3].** For every (y,α) in U, including α = α*, e^{u/z}𝒜 is annihilated by L_pf.

*Proof.* Write ε := e^{𝒲/z}.
1. zθ_yε = −(P/3)ε and z∂_jε = (q_j−α_j+P/3)ε. So for h independent of x_j and of y, q_jhε = z∂_j(hε) + (zθ+α_j)(hε).
2. θ commutes with ∂_x. Iterating step 1 gives Π_j(zθ+α_j)ε ≡ q0q1q2ε modulo total derivatives.
3. θP^m = −(m/3)P^m, so (−3zθ−mz)(P^mε) = P^{m+1}ε. Hence P³ε = −z³Π_j(3θ+j)ε.
4. Since q0q1q2 = yP³, step 3 gives q0q1q2ε = −z³yΠ(3θ+j)ε.
5. Therefore L_pfε dΩ is a sum of total x-derivatives. Apply F2, then F3. ∎

## 4. WKB transport and uniqueness

Let 𝒫 := e^{−u/z}L_pf e^{u/z} = Σ_{m=0}^{3}z^m𝒫_m, with θu = Y.

**Lemma 4 [P].** Let Ω be a connected domain on which Y is holomorphic, L̂(Y) = 0 and YW(Y) ≠ 0. Then:
- (i) 𝒫_0 = 0 and 𝒫_1 = L̂′(Y)θ + ½L̂″(Y)θY + 27yY² = −(W/Y)(θ + ½θ log W). Hence ker 𝒫_1 = C·W^{−1/2}.
- (ii) If f = Σf_kz^k and g = Σg_kz^k have coefficients in O(Ω), 𝒫f = 𝒫g = 0 and g_0 ≢ 0, then f = c(z)g for a unique c ∈ C[[z]]. When α varies, c depends on α only, not on y.
- (iii) Suppose moreover that f and g extend analytically to a boundary point v0 of Ω, with g(v0) ∈ C[[z]] invertible. Then c = f(v0)/g(v0).

*Proof.*
- (i) From (Y+zθ)^nF = Y^nF + z[nY^{n−1}θF + C(n,2)Y^{n−2}(θY)F] + O(z²). Then use L̂′ = −W/Y and θY = −YΠ(Y+α)/W.
- (ii) The z^{k+1} coefficient reads 𝒫_1f_k + 𝒫_2f_{k−1} + 𝒫_3f_{k−2} = 0. 𝒫 commutes with multiplication by series in z alone. If a solution h has h_0 = … = h_{j−1} = 0, then 𝒫_1h_j = 0, so h_j ∈ Cg_0. Subtract c_jz^jg and induct.
- (iii) Evaluate f = c(z)g at v0. ∎

## 5. Normalisation transport

**Proposition 5 [P mod F1–F3, C1, C2].** Fix generic α. Continue R^α_{0i} along any path in the semisimple locus to a point (y,α) and a root Y. Then R^α_{0i} = 𝒜(y,α,Y)/(c·N(z)), where c ∈ C^× is a universal constant fixed by the Gaussian square-root convention. In particular N depends neither on α nor on i.

*Proof.*
1. **Comparison.** By (C1), Lemma 3 and Lemma 4, near y = 0 we have 𝒜_i = κ_i(z)R^α_{0i}, and κ_i does not depend on y.
2. **Coordinates for y→0.** At colour i, q_i = O(y). Put x_i = log y + log p, a translation that adds the constant −α_i log y to u. Then 𝒲_y = yp − α_i log p + Σ_{j≠i}(q_j−α_jx_j) + (pq_jq_k)^{1/3}, which is analytic at y = 0.
3. **Decoupling.** Set w := (pq_jq_k)^{1/3}. This is a linear change in log coordinates, and dΩ = dlog w dx_j dx_k. At y = 0, 𝒲_0 = Σ_{j≠i}(q_j−(α_j−α_i)x_j) + (w − 3α_i log w), which decouples into three nondegenerate one-dimensional critical points.
4. **Limit of 𝒜_i.** By Lemma 2 at y = 0 and Stirling (exmpl:stirling, 4260), lim_{y→0}𝒜_i = c·Π_{a∈{3α_i, α_j−α_i, α_k−α_i}} a^{−1/2}Θ(−z/a).
5. **Limit of R^α_{0i}.** By (C2), and since Δ_i(0) = Π_5 w_i^{(k)}, R^α_{0i}(0) = Π_5 w^{−1/2}Θ(−z/w).
6. **Ratio.** Using Θ(−z/(−ε)) = 1/Θ(−z/ε), κ_i = c·√(ε1ε2)Θ(−z/ε1)Θ(−z/ε2) = cN.
7. **Continuation.** Along any path both sides remain WKB solutions with the same Y, so Lemma 4 on a chain of discs keeps the ratio equal to cN. Because N does not depend on i, any relabelling of colours by monodromy is harmless. ∎

This derives eq:normNLR rather than citing it, and it uses no result of BS lines 4337–4343.

## 6. Continuity to T_F

**Theorem 6 [P mod F1–F3, C1–C3].** Let y0 ∉ {0, −1/27} and let Y* be a root at (y0, α*).
- (a) On a simply connected neighbourhood U of (y0, α*), every coefficient R^{α,[n]}_{0i} extends analytically across the diagonal. Its value there is R^{T_F,[n]}_0 := [z^n]𝒜(·,α*)/(cN). Convergence is locally uniform with all y-derivatives.
- (b) e^{u/z}R^{T_F}_0 solves eq:pf on T_F, and its z^0 coefficient is Δ_i^{−1/2}, with the continued branch.
- (c) The rows R_1, R_2 of eq:RHb, flatness, and Σ_{a,b}R_{ai}(z)η^{ab}R_{bj}(−z) = δ_ij all pass to α = α*. This holds wherever A, D ≠ 0, ∞, so in particular not at y = −1/27.
- (d) For colours 1 and 2, (a) and (b) also hold near y = −1/27.

*Proof.*
- (a) Combine Proposition 5 with Lemma 2.
- (b) Apply Lemma 3 at α*. The z^0 statement follows by continuity of Ψ_{0i} = Δ_i^{−1/2} (C2), which also fixes c.
- (c) The coefficients of eq:RHb do not depend on α except through e1(α), which is analytic. θY = −YΠ(Y+α)/W is analytic on U. The identities hold for generic α (C3) and are analytic.
- (d) The root stays finite and simple, so Lemmas 1–4 apply. ∎

**Definition 7.** The theorem rows are:
- P0^{(i)} := R^{T_F}_{0i};
- P1^{(i)} := C1i(Y_i+zθ)P0^{(i)};
- P2^{(i)} := (Δ_c/C1i²)[(Y_i+zθ)P1^{(i)} − d2P1^{(i)}], with d2 = E(C1i²−1)/(2Δ_c).

Here C1i := 1/c1, θC1i = −C1i·XX and XX := θ log c1. The conifold choice is c1 = C1cf. Frame B multiplies all three rows by h_w(z)^{−1}.

With c1 = A these are exactly eq:RHb, through the dictionary c2 = C1i²/Δ_c = X/A² = D and d2 = e1(α*)E(y). With c1 = A_cf = C1cf, and Cf defined as the identity on G_{T_F} (a convention, made well defined by C4), they are Cf(R_{a i}). Note that buildRmatrix only uses XX = θ log c1. L2b is needed only for Cf to commute with θ on the whole ring.

## 7. Closed form

**Proposition 8 [C: BS 5216–5273, rechecked; P for the constants].** At every colour,
P0^{(i)} = Δ_i^{−1/2}h_w(z)𝔖(v_i,z), and P0^{B,(i)} = Δ_i^{−1/2}𝔖(v_i,z).

*Proof.*
1. **Colour 0, choice of chart.** The critical point is q^{1/3} = c(1,ω,ω²). This is a vector of q^{1/3}-coordinates; in q it is c³(1,1,1). Its entries are distinct, so s_n = e_n(q^{1/3}) is a local biholomorphism there (Disc = −27s3² ≠ 0), and F1 applies.
2. **Colour 0, the computation.** In the s-chart the Hessian is block diagonal and 𝒲 splits exactly. The (s1,s2)-block yields Σr_l(z/s3)^l; this is the definition of r_l. The s3-integral yields (s3^cr)^{−l−1}Θ(−z/E)/Π_{j≤l+1}(1+jz/E), with v = 1/s3^cr. Hence 𝒜_0 = K_0Θ(−z/E)·v𝔖(v,z), where K_0 does not depend on z or y.
3. **Colour 0, constants.** Since N^{−1}Θ(−z/E) = h_w/√(ε1ε2), we get P0^{(0)} = K′h_w v𝔖 with K′ z-independent. The z^0 term must equal Δ_0^{−1/2} = v/√(3ε1ε2E) (Theorem 6(b)), so K′ = 1/√(3ε1ε2E).
4. **Colours 1 and 2.** t ↦ ωt amounts to x_0 ↦ x_0 − 2πi. This shifts 𝒲 by a constant and leaves dΩ unchanged, so the colour-i row is the colour-0 row with t replaced by ω^it. ∎

[E] The following finite-order checks agree with Proposition 8:
- the stored cPP_z9 agrees with the closed form exactly through z^9 at 7 weight pairs;
- z^0 and z^1 agree symbolically (by hand), including h_1 = (E²−ε1ε2)/(12Eε1ε2) and ω_1 = −3/(4E) + 1/(9Et);
- the closed form satisfies q³𝒪f = 0 exactly through z^19;
- the stored rows satisfy eq:pf exactly through z^10 at 6 points.

These checks are not proof steps. The claim that genRmat computes the rows to all orders depends on an unproved ansatz sufficiency and on heuristic reconstruction, and is not used.

## 8. Consequences

**Corollary 9 [P, from Proposition 8].**
- (F2) Δ_0^{1/2}P0^{B,(0)}|_{v=0} = E/(E+z) to all orders. Within the class of Lemma 4(iii) this fixes c(z) ≡ 1.
- (L1) Write ψ = −z and λ_j := (1−jψ/E)^{−1}.
  - [t_c^k]P0^{B,(0)} ∈ span{λ_1,…,λ_k}.
  - [t_c^k](ρ𝔖) ∈ span{λ_1,…,λ_{k+1}}, where ρ = −Ev/t_c.
  - Neither has a constant or a polynomial part.
- (Form) [z^k](v h_w𝔖) is a polynomial in v of degree ≤ k+1. Equivalently, it has the stored shape N_k(Y)/(E+3Y)^{k+1} with deg N_k ≤ k.
- (R) For i = 1, 2, every z-coefficient of P0^{(i)}, P1^{(i)}, P2^{(i)} with c1 = C1cf is analytic at y = −1/27, in either frame.
- (Z∞) lim_{z→∞} z v𝔖(v,z) = −t_c, coefficientwise in v. Equivalently, lim zρ𝔖 = E.

*Proof.*
- **(F2)** Only l = 0 survives at v = 0.
- **(L1)** The l-th summand carries v^{l+1} and the proper fraction z^l/Π_{j≤l+1}(1+jz/E). Its partial fractions lie in span{λ_1..λ_{l+1}}. Since v = t_c·unit with a z-free unit, [t_c^k] receives only l ≤ k−1, or l ≤ k after division by t_c.
- **(Form)** Immediate from the definition of 𝔖.
- **(R)** The coefficients are rational in t and regular where t(1+3ω^it) ≠ 0, and 1−ω^i ≠ 0. C1cf and XX are analytic at Δ_c = 0 with C1cf(0) = −1, and Δ_c·d2 = E(C1i²−1)/2.
- **(Z∞)**
  1. Coefficientwise in v, f := v𝔖 = S_∞/z + O(z^{−2}), with S_∞ = Σr_l(Ev)^{l+1}/(l+1)!.
  2. Write the conjugated operator as 𝒪 = zO_1 + z²O_2 + z³O_3. Its top part is O_3 = Δ_cθ³ + 27yθ² + 6yθ = L2b∘θ.
  3. 𝒪f = 0 holds coefficientwise in v, as identities of rational functions in z. Taking z^{−2}𝒪f as z→∞ gives L2b(θS_∞) = 0.
  4. The exponents at Δ_c = 0 are (0,0). So θS_∞ = κC1cf = κθt_c, and S_∞ = κt_c because S_∞(0) = 0.
  5. The v^1 coefficients (Ev against −Ev) give κ = −1. ∎

**Proposition 10 (symplecticity) [P mod Theorem 6].** For any c1 and XX = θ log c1, and all i, j,
Σ_{a,b}P_a^{(i)}(z)η^{ab}P_b^{(j)}(−z) = δ_ij.

The pairing η^{ab} is read from lp2_lib.m etainvF. Its nonzero entries are η^{00}, η^{11} = −E/k and η^{12} = η^{21} = −1/k, with k = 1/(3ε1ε2) and η^{22} = 0.

*Proof.*
1. **Drop-out identity.** Put Q1 := (Y+zθ)P0 and T := Δ_c(Y+zθ)Q1 + (E/2)Q1. Then P1 = C1iQ1 and P2 = C1i^{−1}[T − zΔ_c·XX·Q1] − C1i(E/2)Q1.
2. In the pairing, the C1i² terms from η^{11} and η^{12} cancel, and so do the XX terms. This leaves Π_ij = η^{00}P0^{(i)}P0^{(j)′} − 3ε1ε2(Q1^{(i)}T^{(j)′} + T^{(i)}Q1^{(j)′}), where ′ means z → −z. So Π_ij depends on P0 alone.
3. At c1 = A, Theorem 6(c) and BS (Symplecticity) give Π_ij = δ_ij. By step 2 the same holds for every c1.
4. **Independent check of the diagonal.** At v → 0 on colour 0: P0 = O(v), Q1(0) = h_w/√(3ε1ε2E), and T(0) = −(E/2)Q1(0). Hence Π_00 = h_w(z)h_w(−z) = 1. This check assumes Π_ii is y-independent, which is round-1 [C]. ∎

**Remark (signs).** The graph sum does not change when the sign of a normalised idempotent is flipped, that is, under (P_{·i}, Δ_i^{1/2}) → (−,−). In every term each occurrence of the index i flips once, and Kronecker terms flip at both ends. The code transports sqDsym together with the rows under t → ω^it, so its branches are compatible.

**Proposition 11 (Givental–Teleman on T_F) [P mod C5].** For g ≥ 2, GW_g^{T_F} equals the graph sum built from R^{T_F}, both at y = 0 and, after continuation, at every y ∉ {0, −1/27}.

*Proof.*
1. Along eq:adelta write GW_g(y,δ) = Σφ_{k,m}(y,δ)F^kA^m (C5). Each φ_{k,m} is a symmetric rational function of the roots, so φ_{k,m} = P_{k,m}/Q with coefficients in C{δ}[y].
2. If δ divided Q, some φ_{k,m}(y0,·) would blow up. This contradicts Theorem 6, under which each graph term is analytic at (y0, α*). So Q(y,0) ≢ 0.
3. Both sides of Q·GW_g = ΣP_{k,m}F^kA^m are regular at δ = 0 (the left side by BS line ~4056), so the identity holds at δ = 0.
4. Pointwise, the limit of the graph sum is the T_F graph sum (Theorem 6). Equality at the remaining points follows by continuity. ∎

**Status of the inputs.**
- Every property the regularity theorem uses is proved here, modulo F1–F3 [S] and the BS genus-0 inputs C1–C3, which the user allows: PF at every colour, uniqueness, the closed form, F2, L1, the stored form, (R), Z∞ and symplecticity.
- The identification with the GW R-matrix in Proposition 11 additionally needs C5.
- F1–F3 are the only items not written out.

The round-2 scripts and their outputs are in `<user home>/AppData/Local/Temp/claude/G--My-Drive-Projects-Old-refGW/37f24460-3601-4c85-b123-5156e19eb475/scratchpad/lp2_round2/` (folders H1_bs, H1_code, continuity, referee_code).