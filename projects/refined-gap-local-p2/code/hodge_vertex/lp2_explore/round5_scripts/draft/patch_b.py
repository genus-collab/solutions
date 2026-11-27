# Round-5 patch, part B: s.2.
PB = []
PB.append((200, 200, 'The closed form of Definition 1.3', r'''The closed form of Definition 1.3 is an algebraic object: all of s.2.1 (PF, Theorem A) and all of s.2.6, **including symplecticity** (Proposition 2.10), are statements about it alone. s.2.2-2.5 and s.2.7 identify it with BS's rows. They are not used in the proofs of Theorems LC_P2 and R_P2, and they are needed for Proposition 2.11 and Proposition 13.2_P2.'''))

PB.append((204, 218, 'Put a := 1 + E/z', r'''Put a := 1 + E/z, e_n := V^n/(a)_n (n in Z, V = E v_i), with (a)_n the Pochhammer symbol. Then v S = z^{-1} Sum_{n>=1} r_{n-1} e_n. L_pf := prod_j (zθ + α_j) + z^3 y prod_{j=0}^{2} (3θ + j) (eq:pf, label at BS 2628, displayed about 2627-2632), and 𝒫 := e^{-u_i/z} L_pf e^{u_i/z} with θ u_i = Y_i (the operator of s.2.4 at T_F).

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
- Finally, Delta_i^{-1/2} = v_i/sqrt(3 ε1 ε2 E), so Delta_i^{-1/2} S is a constant times v S = z^{-1} Sum b_n e_n, and the z-only factors h_w(z), z^{-1} commute with 𝒫. (R3 said "Delta_i^{-1/2} commutes with 𝒫"; it does not, as it depends on y. It is absorbed in v S, as here.) ∎'''))

PB.append((220, 220, '**Corollary 2.2', r'''**Corollary 2.2 (R3:Cor A5, uniqueness).** [P] If b_n in C((z)), b_n = 0 for n <= 0 and 𝒫(Sum b_n e_n) = 0, then b_n = b_1(z) r_{n-1}, so Sum b_n e_n = z b_1(z) v S. Thus PF fixes the row only up to a factor depending on z alone. (Round-3 fix: the factor is z b_1(z), not b_1(z).)'''))

PB.append((222, 222, '### 2.2 Formal', r'''### 2.2 Formal stationary phase (SP1-SP5)

Not used in the proofs of Theorems LC_P2, R_P2. The items were called F1-F5 in R3 and earlier rounds; they are renamed here because F1-F3 are also parts of Theorem LC_P2.'''))

PB.append((251, 251, '**Lemma 2.4', r'''**Lemma 2.4 (R2:Lemma 3, with step 2 written out; R3:Cor B step (i)).** [P, using SP2, SP3] For every (y, α) in U, including α = α*, e^{u/z} 𝒜 is annihilated by L_pf.'''))

PB.append((297, 297, '**Proposition 2.6', r'''**Proposition 2.6 (R2:Prop 5).** [P mod (C1), (C2); uses SP1-SP5] Fix generic α, and let R^α_{0i} be BS's canonical row (eq:canR; (C1): e^{u_i/z} R^α_{0i} formally solves eq:pf; the citation is eq:asymp with eq:canR, not (Flatness)). Continue R^α_{0i} along any path in the semisimple locus (a germ along a path: R^α is defined near y = 0, and the monodromy around y* = u/(27(1-u)) permutes colours on the δ-family, so the statement is about germs, not single-valued holomorphy). Then R^α_{0i} = 𝒜(y, α, Y)/(c N(z)), where c in C^× is a universal constant with c^2 = ±1, fixed by the Gaussian square-root convention. In particular N depends neither on α nor on i.'''))

PB.append((310, 310, '**Theorem 2.7', r'''**Theorem 2.7 (R2:Theorem 6, Prop 8, R3:Cor B; continuity and identification).** [P mod (C1)-(C3); uses SP1-SP5] Let y0 not be in {0, -1/27} and let Y* be a root at (y0, α*).'''))

PB.append((322, 322, '3. **The Wick numbers', r'''3. **The Wick numbers satisfy (Rec) (R3:Cor B).** By Lemma 2.4, e^{u/z} 𝒜_0 solves eq:pf. Θ(-z/E) is invertible and z-only, so 𝒫(v S_Wick) = 0, and Corollary 2.2 with r^W_0 = 1 gives r^W_l = r_l for all l.'''))

PB.append((326, 326, 'Status of step 3', r'''Status of step 3: this makes the closed form S with r_l from (Rec) the same object as BS's Wick expansion. The dependence chain is SP1-SP5 (written) + Lemma 2.4 (written above) + (C1)-(C3).

Weights. Proposition 2.6 and Theorem 2.7 work over BS's generic weights. They are interface statements, and no specialisation of ε from them is used: Proposition 2.10 and everything in s.3-s.6 are proved directly for every complex ε with ε1 ε2 E != 0.'''))

PB.append((333, 335, '**Corollary 2.8', r'''**Corollary 2.8 (Lemma CL, R3 s.0, corrected; carrying identities to the conifold).** [P]
- At every colour, with c1 = C1cf, each z^n coefficient of P^_a^{(i)} is a **Laurent** polynomial in v_i (Lemma 4.1 has 1/v and 1/v^2) multiplied by z-free functions analytic at the conifold (C1i, θv, Δ_c, d2). The same holds in the large-radius frame (the factor h_w depends on z only), and for P_a^{B,(i)} = Delta_i^{-1/2} P^_a^{(i)} up to the constant sqrt(3 ε1 ε2 E)^{-1}. Note Delta_0^{±1/2} is proportional to v_0^{∓1}, so it is a factor of this Laurent type and is not analytic. So each such coefficient is a meromorphic germ at y = -1/27, and an identity between such series that holds on a punctured neighbourhood holds as an identity of Laurent expansions in t_c.
- **At colour 0 only, and in frame B only:** since v_0 = O(t_c), each [t_c^k] coefficient of ρ^e P^^{(0)} involves only finitely many v-orders. So [z^n] and [t_c^k] commute, and each [t_c^k] coefficient is a rational function of z. This fails at colours 1, 2, where v_j(0) != 0 and [t_c^k] of S(v_j, z) is only a formal power series in z, and it fails in the large-radius frame, since h_w is not rational. s.3.1 and Proposition H2 use only the colour-0, frame-B case.'''))

PB.append((338, 338, '- (F2) S(0, z)', r'''- (F2) S(0, z) = E/(E + z): only l = 0 survives at v = 0.'''))

PB.append((348, 348, '*Proof of (Z-infinity).*', r'''*Proof of (Z-infinity).* Coefficientwise in v, f := v S = S_∞/z + O(z^{-2}), with S_∞ = Sum r_l (Ev)^{l+1}/(l+1)!. Write the conjugated operator as 𝒪 = zO_1 + z^2 O_2 + z^3 O_3, with top part O_3 = Δ_c θ^3 + 27 y θ^2 + 6 y θ = L2b∘θ. The V-adic form of Theorem A gives 𝒪f = 0 coefficientwise in v, as identities of rational functions in z. Taking z^{-2} 𝒪 f as z -> ∞ gives L2b(θ S_∞) = 0. The exponents at Δ_c = 0 are (0,0), so θ S_∞ = κ C1cf = κ θ t_c, and S_∞ = κ t_c because S_∞(0) = 0. The v^1 coefficients (E v against -E v, Lemma 1.2(b)) give κ = -1. ∎'''))

PB.append((352, 363, '**Proposition 2.10', r'''**Proposition 2.10 (symplecticity, K = 1, orthogonality; algebraic proof, round 5).** [P] Let ε1, ε2 be complex numbers with ε1 ε2 E != 0, or indeterminates. For any non-vanishing section c1 with XX = θ log c1, and all i, j,

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
- This discharges the item marked [O, optional] in round 4 (an algebraic proof of Π_ij = δ_ij from the closed form). It has had one referee pass (sound with minor fixes, applied here; the referee re-derived Theorem A and S0-S6 by hand).'''))

PB.append((365, 365, '**Proposition 2.11', r'''**Proposition 2.11 (R2:Prop 11; Givental-Teleman on T_F).** [P mod (C1)-(C3), (C5)] Not used in Theorems LC_P2, R_P2. For g >= 2, GW_g^{T_F} equals the graph sum built from R^{T_F} at every y not in {0, -1/27}, by continuation. *Proof.* Along eq:adelta (3232), GW_g(y, δ) = Sum φ_{k,m}(y, δ) F^k A^m (C5), with φ_{k,m} symmetric rational functions of the roots, φ_{k,m} = P_{k,m}/Q with coefficients in C{δ}[y]. If δ divided Q, some φ_{k,m}(y0, ·) would blow up, contradicting Theorem 2.7 (each graph term is analytic at (y0, α*)). Both sides of Q·GW_g = Sum P_{k,m} F^k A^m are regular at δ = 0 (the left by BS line ~4056), so the identity holds at δ = 0; the pointwise limit of the graph sum is the T_F graph sum. ∎ Note: BS line ~4056 is argued only briefly; see s.8.3. BS's prop:fingenR (3611) is stated for rational tori only; at T_F, R_0's independence of A follows from Theorem 2.7(e), since the row is rational in t.

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

**Scope.** Propositions 2.10 and 2.12 reprove BS's (Symplecticity) on T_F independently of eq:sympl and of Theorem 2.7(c), but not independently of BS altogether: identifying row 0 of Definition 1.3 with BS's R_0 still needs Theorem 2.7(e), hence SP1-SP5, (C1) and (C2). The proof of Theorem 2.7(e) uses neither eq:sympl nor Theorem 2.7(c), so there is no circularity. BS lem:modrelorb (about 4581) is also proved only by a sketch; it is outside this chain and was not examined.'''))
