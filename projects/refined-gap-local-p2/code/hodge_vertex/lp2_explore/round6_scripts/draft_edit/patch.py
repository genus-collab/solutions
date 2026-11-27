# Round 6 targeted edits to LP2_REGULARITY_PROOF_draft.md. Every 'old' must occur exactly once.
import os
import sys
F = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "LP2_REGULARITY_PROOF_draft.md")
s = open(F, encoding="utf-8", newline="").read()
fails = []

def rep(old, new):
    global s
    c = s.count(old)
    if c != 1:
        fails.append((c, old[:90]))
        return
    s = s.replace(old, new)

def span(start, end, new, keep_end=True):
    """Replace from the unique marker start up to the first following end marker."""
    global s
    c = s.count(start)
    if c != 1:
        fails.append((c, start[:90]))
        return
    i = s.index(start)
    j = s.find(end, i)
    if j < 0:
        fails.append((0, end[:90]))
        return
    k = j if keep_end else j + len(end)
    s = s[:i] + new + s[k:]

# ---- header
rep('round 5 wf_9d7fd94e-8d7. NOT yet refereed as a whole. -->',
    'round 5 wf_9d7fd94e-8d7. Round 6 (2026-11-02): edited in place from LP2_REGULARITY_PROOF_draft_round5.md, applying the referee reports D_rows (s.0-2), D_core (s.3-6), D_hodge (s.7.0-7.7), D_VI (s.7.8) and the items of L_logic (LEMMA_A_PROOF.md) that touch this file; all "sound with minor fixes", no major item. Not refereed as a whole: see s.0.4. -->')

# ---- s.0.1 colour remark (D_rows)
rep('it is the unique colour whose critical point escapes (Lemma 2.1).',
    'it is the unique colour with Y_i -> ∞ (the critical point escapes), since Δ_c = (1+3t)(1-3t+9t^2), and at t = -1/3 one has 1 + 3ω^i t = 1 - ω^i != 0 for i = 1, 2.')

# ---- s.0.2 labels
rep('**Theorem R_P2 (regularity of graphs with an edge).** [P mod Lemma A(h), 2 <= h <= g-1]',
    '**Theorem R_P2 (regularity of graphs with an edge).** [P mod (MV), through Lemma A_P2(b) for 2 <= h <= g-1; (MV) is the Marino-Vafa formula, s.0.3]')
rep('The same therefore holds for the colour sum.\n',
    'The same therefore holds for the colour sum.\n\n**M_irr is discharged.** M_irr(h) holds for every h >= 2 by Lemma A_P2(b) (s.7.4; P modulo (MV) and textbook facts). So the hypothesis of Theorem R_P2 and of Corollary SV_P2 is not an open assumption. It is kept in their statements only to show which h are used.\n')
rep('[P mod Lemma A(h) for 2 <= h <= g, and Theorem VI, which is P mod LLZ eqs. (1), (5) and the boson-fermion correspondence (s.7.8)]',
    '[P mod (MV), through Lemma A_P2(b) for 2 <= h <= g, and through Theorem VI, which needs LLZ eq. (1) only, and only at μ = 1^a, plus textbook facts (string, dilaton, Frobenius, Murnaghan-Nakayama, principal specialisation) and (7.2) (s.7.8)]')

# ---- s.0.3
rep('(C4) prop:algindABC 2854, only to make Cf well defined; and t_cf',
    '(C4) prop:algindABC 2854, to make Cf well defined and for the step of Proposition 2.11 from "δ divides Q" to "some φ_{k,m} blows up" (a regular graph sum does not bound the individual φ_{k,m} without independence of A, F); and t_cf')
rep('proved there modulo the Marino-Vafa formula (LLZ eq. (1); OP Theorem 1 is a second published proof) and textbook facts. The assembled file has not yet had its own referee; see s.7.4 and s.8.2.',
    'proved there modulo the Marino-Vafa formula, called (MV) below (LLZ eq. (1); OP Theorem 1 is a second published proof) and textbook facts. The assembled file had its own referee pass in round 6 (L_logic: sound with minor fixes, no major item); see s.7.4 and s.8.2.')
rep('from LLZ eqs. (1) and (5), the boson-fermion correspondence, the string and dilaton equations and (7.2). It no longer',
    'from LLZ eq. (1), only at μ = 1^a, textbook facts (string, dilaton, Frobenius, Murnaghan-Nakayama, principal specialisation) and (7.2). It no longer')

# ---- s.0.4
rep('- That the whole text has been refereed. Proposition 2.10, Proposition 2.12 and Theorem VI have had one referee pass each (sound with minor fixes, applied here); the extension of Fpot_h in s.1.5 and (H0) of s.7.1 are new in this edit and unrefereed.',
    '- That the whole text has been refereed. Round 6 gave each assembled part one referee pass of its own, each with verdict sound with minor fixes and no major item; the fixes are applied here. Coverage: s.0-s.2 apart from the items below, including the extension of Fpot_h in s.1.5 and (H0) of s.7.1 (D_rows); s.3-s.6 (D_core); s.7.0-s.7.7, with (H0) re-derived again (D_hodge); s.7.8 (D_VI); LEMMA_A_PROOF.md (L_logic). Proposition 2.10, Proposition 2.12 and Theorem VI also had a pass in round 5. **Not refereed:** s.2.2-2.5 including Theorem 2.7, Proposition 2.11, Remark 2.13(a)-(c), BS\'s unproved assertion behind Lemma 1.2(c) (float check only; used only for Corollary G), Lemma 1.4 at genus 4, and s.8.3. D_rows and D_hodge judged by reading that s.2.2-2.5 and s.2.7 are not used in the proofs of Theorems LC_P2, R_P2 and Corollary SV_P2. All finite-order computations of those passes are evidence only.')

# ---- s.1.1 notation (D_rows)
rep('are written e_k(α), never e_k alone.',
    'are written e_k(α), never e_k alone (exceptions, flagged where they occur: the functions e_n := V^n/(𝖺)_n of s.2.1, and the abbreviation e_k := e_k(α) in Remark 2.13).')

# ---- Lemma 1.2(c), remark (D_rows, D_core)
rep('Vanishing kills the constant, and the linear coefficient gives κ = 1. ∎',
    'Vanishing kills the constant, and the linear coefficient gives κ = 1. On y in (-1/27, 0) the expression at BS 5145 is meant with log y = log|y| + iπ; with the other branch it does not vanish at Δ_c = 0 [E, float at 5 values of y]. ∎')
rep('and agree with the proof, but it has had no separate referee pass. It is used only for Corollary G.',
    'and agree with the proof. BS asserts the value of the combination at 5145 without proof; round 6 (D_rows) tested it only by float [E]. It is used only for Corollary G.')
rep('All filtration statements hold for any ρ(0) != 0.',
    'The filtration statements (the memberships of L1-L3) are unchanged if t_c is rescaled as an expansion variable, t_c -> a t_c with c1 = C1cf kept, which multiplies ρ by the constant 1/a; then ρ(0) = 1/a and the base point is q_0/a. 𝔇 and (H0) as written, hence s.7.5(i)-(ii), need ρ(0) = 1. This is not a change of c1: c1 = a θ t_c gives C1i(0)^2 != 1, and d2 acquires a pole. The statements are also not stable under a non-linear change of coordinate with the same leading term: with Δ_c in place of t_c in ρ, the constant of [t_c^1](ρ q_t) is (11/18) E != 0, because Z-infinity is the exact identity S_∞ = -t_c.')

# ---- s.1.3 Cf (D_rows)
rep('With c1 = C1cf = A_cf, and Cf taken **by convention** to be the identity on G_{T_F} (well defined by (C4); BS sec:cfgen describes Cf only through τ -> -1/(3τ), which does not fix the Hauptmodul), they are Cf(R_{ai}).',
    'With c1 = C1cf = A_cf they are Cf(R_{ai}). Here Cf is BS\'s ring homomorphism G_{T_F}[A,F] -> G_{T_F}[A_cf,F_cf] induced by Cf(A) = A_cf, Cf(F) = F_cf (BS lines 5148-5160; well defined by (C4)), so it is the identity on G_{T_F}, and D_cf := X/A_cf^2 treats X as fixed. BS\'s further gloss Cf(f)(τ) = f(-1/(3τ)) cannot be literal on G_{T_F}, since it does not fix the Hauptmodul.')

# ---- s.1.5 (D_rows)
rep('Let 𝐑 := C[[t_c]] (or Q(ε1, ε2)[[t_c]], or any t_c-adically complete Q-algebra of coefficients).',
    'Let 𝐑 := 𝐑_0[[t_c]], with 𝐑_0 a Q-algebra in which ε1, ε2, E are invertible (for instance C or Q(ω, ε1, ε2)).')

# ---- s.2 opening (D_rows)
rep('all of s.2.1 (PF, Theorem A) and all of s.2.6, **including symplecticity** (Proposition 2.10), are statements about it alone.',
    'all of s.2.1 (PF, Theorem A) and all of s.2.6 except Proposition 2.11, **including symplecticity** (Proposition 2.10), are statements about it alone.')

# ---- s.2.1 (D_rows): rename a -> 𝖺, Q(E)
rep('Put a := 1 + E/z, e_n := V^n/(a)_n (n in Z, V = E v_i), with (a)_n the Pochhammer symbol.',
    'Put 𝖺 := 1 + E/z, e_n := V^n/(𝖺)_n (n in Z, V = E v_i), with (𝖺)_n the Pochhammer symbol.')
rep('V e_n = (a+n) e_{n+1} and z(a+n-1) = E + nz give',
    'V e_n = (𝖺+n) e_{n+1} and z(𝖺+n-1) = E + nz give')
rep('Since 1/(a)_n = prod_{j<=n} z/(E + jz), each e_n with n >= 0 lies in V^n z^n Q[[z]], so',
    'Since 1/(𝖺)_n = prod_{j<=n} z/(E + jz), each e_n with n >= 0 lies in V^n z^n Q(E)[[z]], so')
rep('If b_n in C((z)), b_n = 0 for n <= 0 and 𝒫(Sum b_n e_n) = 0, then',
    'If b_n in C((z)), b_n = 0 for n <= 0 and 𝒫(Sum b_n e_n) = 0 (the sum understood V-adically, since b_n in C((z)) need not converge z-adically), then')

# ---- Cor 2.8, 2.9 (D_rows, D_core)
rep('z-free functions analytic at the conifold (C1i, θv, Δ_c, d2).',
    'z-free functions analytic at the conifold (C1i^{±1}, XX, Δ_c, d2).')
rep('They are rational in t and regular where t(1 + 3ω^i t) != 0',
    'They are polynomials in C1i^{±1}, XX, d2 with coefficients rational in t, regular where t(1 + 3ω^i t) != 0')
rep('Taking z^{-2} 𝒪 f as z -> ∞ gives L2b(θ S_∞) = 0.',
    'O_1 and O_2 have poles at V = 0 (e.g. O_1 = -(W/Y)(θ + ½ θ log W), s.2.4), so everything is read coefficientwise in Laurent series in V: each V-coefficient of 𝒪 f is a finite combination of V-coefficients of f, with coefficients rational in z, and each V-coefficient of f is O(1/z). So only z^3 O_3 contributes in the limit, and taking z^{-2} 𝒪 f as z -> ∞ gives L2b(θ S_∞) = 0. Equivalently: e_n -> V^n/n!, 𝗌 -> d/dV, 𝖭 -> V d/dV, 𝗍 -> 3θ, so 𝒫/z^3 -> θ^3 - (d/dV)^3 = L2b∘θ.')

# ---- Prop 2.10 (D_rows)
rep('*Proof.* Work over K := C(t), y = t^3, θ = (t/3) d/dt, with coefficients series in z. Write Th_i',
    '*Proof.* Let 𝔽 be the constant field (C, or Q(ω, ε1, ε2)) and K := 𝔽(t), y = t^3, θ = (t/3) d/dt, with coefficients series in z. S0 takes place in a differential extension of K containing C1i^{±1}, XX and a square root of 3 ε1 ε2 E; S1-S6 take place in K[[z]] only. Write Th_i')
rep('It has had one referee pass (sound with minor fixes, applied here; the referee re-derived Theorem A and S0-S6 by hand).',
    'It has had a referee pass in round 5 (sound with minor fixes, applied here; the referee re-derived Theorem A and S0-S6 by hand) and was re-derived by hand again in round 6 (D_rows, D_core).')

# ---- Prop 2.11 (C4), consequence of the s.0.3 fix
rep('**Proposition 2.11 (R2:Prop 11; Givental-Teleman on T_F).** [P mod (C1)-(C3), (C5)]',
    '**Proposition 2.11 (R2:Prop 11; Givental-Teleman on T_F).** [P mod (C1)-(C5)]')
rep('If δ divided Q, some φ_{k,m}(y0, ·) would blow up,',
    'If δ divided Q, then, by the algebraic independence of A, F over functions rational in y ((C4)), some φ_{k,m}(y0, ·) would blow up,')

# ---- Remark 2.13(c) (D_rows)
rep('eq:XFABC (2857)', 'eq:XFABC (2860)')

# ---- s.3 header, E-scaling (D_core)
rep('Throughout s.3-s.5, v := v_0, on the real branch, c1 = C1cf. Then θv = 1/E - v/3, Y = -E/3 + 1/v, θ log v = Y/E, Delta_0 ∝ v^{-2}.\n',
    'Throughout s.3-s.5, v := v_0, on the real branch, c1 = C1cf, except Lemma 4.1, Lemma 4.2(G), H2(d) and Proposition 3.1, which allow general c1. Then θv = 1/E - v/3, Y = -E/3 + 1/v, θ log v = Y/E, Delta_0 ∝ v^{-2}.\n\n'
    '**Remark (E-scaling).** [P] V = E v, ρ, t_c and κ (Lemma 4.2) do not depend on ε; S(v, z) = G(V, z/E) with G(V, u) := Sum_l r_l (Vu)^l / prod_{j<=l+1} (1 + ju); and Y, Ŷ, μ are E times ε-free series. Hence σ_a(ψ; t_c; E) = E^a σ_a(ψ/E; t_c; 1), likewise for the colour-j rows, and\n\n'
    '    Delta_i = ε1 ε2 E^3 × (universal),   Dil^_i(d) = E^{1-d} × (universal),   Ehat_ij = ε1 ε2 E^2 Uhat_ij(ψ/E, w/E; t_c), with Uhat_ij universal.\n\n'
    'Consequences.\n'
    '- Theorem LC_P2 is one statement over Q (Q(ω) at the regular colours). The framing τ = -ε1/E enters the rows only through the prefactor ε1 ε2 and the scale E; τ matters only on the Hodge side.\n'
    '- Finite-order row and edge checks at several weight pairs are a single test of the scale-free object, plus a test of the scaling (Appendix). They do not probe non-special framing.\n'
    '- Unlike the conifold, whose top leaf vanishes at c = (k+1)/i (LSP Theorem LC), the top coefficient of σ_0^(k) is non-zero iff r_k != 0, independently of ε, since only the summand l = k of S reaches λ_{k+1}. r_k is non-zero with oscillating sign through k = 200 [E]; the pattern starts +++++------++++++---. So the sharpness clause of LSP Theorem LC has no literal analogue here, and none is claimed.\n')

# ---- s.3.1 (D_core)
rep('By Corollary 2.8 (the colour-0, frame-B case), the σ^(k) lie in the space V := C[ψ]',
    'By Corollary 2.8 (the colour-0, frame-B case) the σ^(k) are rational in ψ; by Corollary 2.9(L1) and (I1), (I2) of Lemma 4.2 they lie in the space V := C[ψ]')
rep('as an identity of Laurent expansions in t_c with coefficients in V ⊗ V. This is',
    'as an identity of Laurent expansions in t_c whose coefficients are rational functions of ψ, namely the restrictions to ψ\' = -ψ of elements of V ⊗ V. This is')
rep('no continuation or branch argument is needed; Corollary 2.8 gives the coefficients in V.',
    'no continuation or branch argument is needed.')

# ---- s.4 (D_core)
rep('(I1) σ_1 = D_c σ_0, where Ŷ := C1i Y = E d log v/dt_c = E/t_c + E/18 + O(t_c);',
    '(I1) σ_1 = D_c σ_0, where Ŷ := C1i Y = E d log v/dt_c = E/t_c + E/18 + O(t_c) (Ŷ = E/t_c + E d log ρ/dt_c, by v = -ρ t_c/E);')
rep('and (LAM2) holds. Here Ŷ^(m) := [t_c^m] (Ŷ - E/t_c).',
    'Here Ŷ^(m) := [t_c^m] (Ŷ - E/t_c). *Proof.* With z = -ψ, the coefficient of t_c^n of (E - ψ) F/t_c - ψ dF/dt_c is (E - ψ - (n+1) ψ) F^(n+1). ∎')
rep('With r_1 = 1/9 and v = -t_c/E - t_c^2/(18E) + O(t_c^3): σ_0^(1)',
    'With r_1 = 1/9 and v = -t_c/E - t_c^2/(18E) + O(t_c^3) [from Z-infinity: with V = E v, S_∞ = V + V^2/18 + O(V^3) = -t_c, so ρ = 1 + t_c/18 + O(t_c^2)]: σ_0^(1)')
rep('The map between the two gauges is z-free, analytic and invertible (C1cf(0) = -1), so H2 holds in either gauge. In particular',
    'The map between the two gauges is z-free, analytic and invertible (C1cf(0) = -1): P_1^hy = c1 P_1, P_2^hy = C1i P_2 + (E/2)(C1i^2 - 1) c1 P_1. It preserves the membership in C + Λ_m, but not H2(a), (c) or the no-constant clause of H2(b): a constant at k >= 0 appears, e.g. -(2/9)E in σ_1^hy at k = 0, and σ^hy(-1) = (0, -E, E^2). This is harmless, since Ehat_ij kills every ψ-constant and H2 is used only in the true gauge. In particular')

# ---- s.5 (D_core)
rep('**Theorem E00 (R3:Theorem E00).** [P, using s.3.1 and H2 Steps 1-2] By (1.1) and K = 1, Ehat_00 = [Delta_0 - ⟨τ(ψ), τ(ψ\')⟩]/(ψ + ψ\'), τ(ψ) := P^(-ψ).',
    '**Theorem E00 (R3:Theorem E00).** [P, using s.3.1 and H2 Steps 1-2; (iv) also uses σ_0^(0), σ_1^(0) of H2 Step 4] By (1.1) and K = 1, Ehat_00 = [Delta_0 - ⟨𝗉(ψ), 𝗉(ψ\')⟩]/(ψ + ψ\'), 𝗉(ψ) := P^(-ψ).')
rep('At n = 0, 1 there are no admissible pairs, giving (ii).',
    'At n = 0, 1 there are no admissible pairs. This gives the formal statement of (ii). By Corollary 2.8 each Ehat_00(m,n), a polynomial in finitely many z-coefficients of the rows, is a meromorphic germ at t_c = 0, so having no negative powers it is analytic.')

# ---- s.7.0 (D_rows, D_core)
rep('colour 0 (singular on the real branch; Lemma 2.1)',
    'colour 0 (singular on the real branch: 1 + 3ω^i t -> 1 - ω^i != 0 for i = 1, 2)')
rep('    ρ = τ/(t alg)                          ρ = -E v_0/t_c\n',
    '    ρ = τ_LSP/(t alg), τ_LSP := 1 - e^-t   ρ = -E v_0/t_c\n')
rep('    w = (c, -x, -1)                        w = (E, -ε1, -ε2) = ε2 (c, -x, -1) at x = ε1/ε2, using c = 1 + x (LSP App. A)\n',
    '    w = (c, -x, -1)                        w = (E, -ε1, -ε2) = ε2 (c, -x, -1) at x = ε1/ε2, using c = 1 + x (LSP App. A)\n'
    '    (LAM2), (LAM3), (LAM6)                 (LAM1), (LAM2), (LAM3)\n'
    '    L_i at colour 0                        Ŷ\n'
    '    r_{1,i} = (L_i + z D_i) r_{0,i},       P^_1 = (Ŷ + z d/dt_c + z d log v/dt_c) S, which is (I1) in LSP\'s form\n'
    '      D_i = d/dt - Delta_i\'/(2 Delta_i)\n')

# ---- s.7.2 (D_core)
rep('it never reaches a vertex.\n',
    'it never reaches a vertex.\n\nH2 Step 3 (the constants of σ_2) is not load-bearing for E00 or Theorem M, because constants of σ feed only C_i0, C_0l, C_00 and cancel in (M-a). What is load-bearing is Z-infinity, in two places: H2 Step 1 (no constant in σ_1 at k >= 0, hence no ψ-linear term in σ_2) and Lemma D0(iii).\n')

# ---- s.7.3 (D_hodge)
rep('and dim >= 2 nd, so a genus-0 vertex with m = 0 needs at least 3 half-edges.',
    'and dim >= 2 nd, so a genus-0 vertex with n_v <= 2 half-edges could carry no non-zero decoration for any m (dim = n_v + m - 3 < 2m); this agrees with Γ being stable.')

# ---- s.7.4
rep('[C: hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md, Theorem L and Corollaries L1-L3, s.7.5 there, for all h >= 2; there P modulo (MV) and textbook facts; status below]',
    '[P mod (MV) and textbook facts: hodge_vertex_conifold/lemmaA/LEMMA_A_PROOF.md, Theorem L and Corollaries L1-L3, s.7.5 there, for all h >= 2; status below]')
rep('(eqs. (1)-(3) of arXiv:math/0306434v3; published JDG numbering not confirmed; "LLZ Thm 1"',
    '(eqs. (1)-(3) of LLZ, J. Differential Geom. 65 (2003) 289-340, p. 291, numbered as in arXiv:math/0306434v3; "LLZ Thm 1"')
rep('The assembled LEMMA_A_PROOF.md has not had its own referee (its open point O1). That file\'s status table',
    'The assembled LEMMA_A_PROOF.md had its own referee pass in round 6 (L_logic: sound with minor fixes, no major item), which answers its open point O1. That file\'s status table')
rep(' and the empty profile, which is Theorem VI here. Its open point O5 and its s.8.5, written in parallel with this edit, still list (VI) and lem:modrel as resting on R2 Cor. 7 and on critic 1 S4; s.7.8 and Proposition 2.12 of this text supersede that.',
    ' and the empty profile, which is Theorem VI here.')

# ---- s.7.5 (D_hodge)
rep('Dilaton decorations at every vertex are summed first, as in (1.2).',
    'Dilaton decorations at every vertex are summed first, as in (1.2) (for g_v <= 1 this is just the finite sum of (0.1); (1.2) and Fpot are defined for h >= 2 only).')
rep('and they are read coefficientwise in the edge variables of the far end.\n',
    'and they are read coefficientwise in the edge variables of the far end. Precisely: 𝔡^s is 𝐑-multilinear and depends only on a finite jet of each slot (m <= 3h-3+s, and the slot degree is cut off by the dimension), so it is t_c-adically continuous in each slot and may be expanded term by term over t_c-series of leaves. At an unrescaled far end the exponent b is bounded by 3g_u-3+n_u (s.0.1), so only finitely many [w^b] occur. At fixed t_c-order, the number r of δ-insertions is at most that order, and only finitely many leaf multisets R occur.\n')
rep('Only two properties of ρ are used: it is an analytic z-free unit, and ρ^2 Delta_0 is a constant times t_c^{-2}.',
    'Three properties of ρ are used: it is analytic and z-free; ρ(0) = 1, so that λ = ρ is admissible in (H0) (λ in 1 + t_c 𝐑) and [t_c^0](ρ q_t) = q_0, the base point of 𝔇; and ρ^2 Delta_0 = 3 ε1 ε2 E^3 t_c^{-2}.')
rep('- their dilaton and edge factors are analytic (Lemma D0(iv), (M-e), E00(ii));\n- the unit ρ^{-1} does not raise leaf degree.\n',
    '- their dilaton and edge factors are analytic (Lemma D0(iv), (M-a), (M-d), (M-e), E00(ii)).\n')
rep('All leaves at order t_c^k are λ_j with 1 <= j <= k + 1: only lattice points of M enter.',
    'All leaves at order t_c^k are λ_j with 1 <= j <= k + 1.')
rep('- If R\' is non-empty, D(R\') = D(R) <= k. If k <= 2h - 3 the term vanishes by M_irr(h), with h <= g - 1 by Lemma 7.1.',
    '- If R\' is non-empty, then D(R\') = D(R), and the term vanishes by M_irr(h) whenever D(R) <= 2h - 3, whatever the t_c-order k at which it occurs (h <= g - 1 by Lemma 7.1).')
rep('per colouring, for every genus, modulo Lemma A_P2(b) for 2 <= h <= g-1. Genus 2',
    'per colouring, for every genus, modulo Lemma A_P2(b) for 2 <= h <= g-1, i.e. modulo (MV). Genus 2')

# ---- s.7.6
rep('**Theorem VI (vertex identity (VI)).** [P mod LLZ eqs. (1), (5), the boson-fermion correspondence, string/dilaton; proof in s.7.8]',
    '**Theorem VI (vertex identity (VI)).** [P mod LLZ eq. (1), only at μ = 1^a, plus textbook facts (string, dilaton, Frobenius, Murnaghan-Nakayama, principal specialisation) and (7.2); proof in s.7.8]')
rep('= (2g-3)! P_g(x,1)/x. This is conifold F4 at n = 0.',
    '= (2g-3)! P_g(x,1)/x, with S(z) := (z/2)/sinh(z/2) and P_g(a,b) := [t^{2g}] S(at) S(bt). This is conifold F4 at n = 0.')
rep('The proof in s.7.8 uses no TR, no Eynard/ABDKS, no Lemma A, no Teleman and no finite-order input.',
    'The proof in s.7.8 uses no TR, no Eynard/ABDKS, no Theorem L (it quotes only LEMMA_A_PROOF.md Proposition 2.1, Lemma 3.2 and eq. (4.2), whose proofs do not involve Theorem L), no Teleman and no finite-order input.')
rep('See also BS prop:cfgraphsum (statement at 5114, eq:betagiv 5122, label 6366), which is related and which has not been reconciled with the single-vertex value Z[g,0].',
    'BS prop:cfgraphsum (statement at 5114, eq:betagiv 5122, label 6366) is (VI) in Givental form. Since Υ(-z/ε1, -z/ε2) = h_w(-z) = h_w(z)^{-1}, the graph sum of BS (5096-5124), with pure ψ vertices, edge [1 - h_w(z)^{-1} h_w(w)^{-1}]/(z+w) and dilaton leaf 1 - h_w(z)^{-1}(1 - z/E)^{-1}, is by Mumford (Omega_w = h_w·1) the single-vertex value Fpot_g(q_0) = Z[g,0], because q_0 = -z(1 - z/E)^{-1}. For general g this reduction is a derivation sketch; it was tested exactly only at g = 2, (ε1, ε2) = (7/3, 5/2), all 7 graphs: BS sum = -152041/4097352000 = Z[2,0], and this times -E^3 ε1 ε2 is β_2 = 152041/6220800 [E]. BS\'s proof runs through prop:ResConifold, so it stays a cross-check and does not replace s.7.8.')

# ---- s.7.7(c) (D_rows)
rep('Near y = -1/27 the roots t_i are analytic and distinct (Lemma 2.1), so labels',
    'Near y = -1/27 the roots t_i = ω^i t are analytic and distinct, since t != 0, so labels')

# ---- s.7.8 Inputs
span('**Inputs.** LLZ eq. (1)-(3) and Theorem 1 eq. (5)',
     ' The proof combines the MV formula,',
     '**Inputs.** LLZ eq. (1)-(3) (J. Differential Geom. 65 (2003) 289-340 and arXiv:math/0306434v3, read), used only at μ = 1^a; LEMMA_A_PROOF.md Proposition 2.1 (conversion to real form), Lemma 3.2 (dilaton point) and s.4 eq. (4.2) (first mKP equation, derived there), none of which involves Theorem L; textbook facts: string and dilaton equations for pulled-back classes, Int_{1,1} ψ = Int_{1,1} λ_1 = 1/24 (which fixes κ_0), Frobenius\' transposition identity (classical; no number pinned), Murnaghan-Nakayama and the Frobenius formula (inside LEMMA_A_PROOF.md s.4), and the principal specialisation s_ν(q^ρ) = q^{κ_ν/4}/prod_box [h] (Macdonald I.3 Ex. 2, recalled, example number not confirmed; LLZ eq. (12) states it; the exponent was re-derived by hand); (7.2).')

# ---- s.7.8 Step 1 (shrink: LEMMA_A Prop 2.1)
span('*Step 1 (the MV formula in real form, all g and all μ).*',
     'By the principal specialisation, 𝖹 = Sum_ν',
     '*Step 1 (the MV formula in real form).* [C LLZ eqs. (1)-(3); the conversion is P: LEMMA_A_PROOF.md Proposition 2.1, eq. (2.1), which writes out the hook form of V_ν, the substitution x = iλ and every sign, uniformly in g and μ] Put 𝖹(p) := Sum_ν e^{(τ+1/2) κ_ν x/2} prod_box (2 sinh(hx/2))^{-1} s_ν(p) and 𝖥 := log 𝖹. With l = l(μ), Π_μ = prod_i Π_{μ_i}, Π_m = prod_{a=1}^{m-1}(mτ + a)/(m-1)!, M_g(μ) := ⟨Omega prod_i 1/(1-μ_i ψ_i)⟩_{g,l} and M_0(μ) := |μ|^{l-3},\n\n'
     '    [x^{2g-2+l} p_μ] 𝖥 = σ^{g+l-1} Π_μ M_g(μ)/|Aut μ|.   (B1)\n\n'
     'Only μ = 1^a is used below: ∂_2 𝖥 is eliminated in Step 2(c) without (B1). ')

# ---- s.7.8 Step 2(a) (D_VI)
span('(a) Cut-and-join. [C: LLZ Theorem 1 eq. (5)',
     '    ∂_τ 𝖥|_p = x p_1^2 ∂_2 𝖥.',
     '(a) Cut-and-join on the slice. [P, Frobenius] For ν |- n, s_ν = Sum_μ χ_ν(μ) p_μ/z_μ gives at p = (p_1, 0, ...): s_ν = χ_ν(1^n) p_1^n/n! and ∂_2 s_ν = χ_ν(2,1^{n-2}) p_1^{n-2}/(2(n-2)!). Frobenius\' classical identity binom(n,2) χ_ν(2,1^{n-2}) = f_2(ν) χ_ν(1^n) gives p_1^2 ∂_2 s_ν = f_2(ν) s_ν on the slice (for n <= 1 both sides vanish). Hence ∂_τ 𝖹|_p = x p_1^2 ∂_2 𝖹 on the slice, and dividing by 𝖹 gives\n\n')

# ---- s.7.8 Step 2(b) (D_VI)
span('(b) First mKP. [P, with C for the boson-fermion dictionary]',
     '    2 ∂_2(𝖥 - 𝖥_1) + ∂_1^2(𝖥 + 𝖥_1)',
     '(b) First mKP. [P: LEMMA_A_PROOF.md s.4] Put t_k = p_k/k, β := τ x, A = e^β, p\'_k = 1/(2 sinh(kx/2)). By LEMMA_A_PROOF.md (F2) and its Specialisation paragraph, Z_0 = 𝖹 and Z_1 = 𝖹(e^{kβ} p_k), with no normalising constant, and (4.2) holds coefficientwise in Q[τ]((x))[[p]]. At k = 2, 2P_2 = H_1^2 + H_2, where Z_0(b+e) Z_1(b-e) = exp(Sum e_k H_k) Z_0 . Z_1. So (H_1^2 + H_2) Z_0 . Z_1 = 0. Both series have constant term 1 in p, so one may divide by Z_0 Z_1. Since ∂/∂t_1 = ∂_1 and ∂/∂t_2 = 2∂_2, and with 𝖥_1(p) := 𝖥(e^{kβ} p_k), this reads:\n\n')
rep('    2 ∂_2(𝖥 - 𝖥_1) + ∂_1^2(𝖥 + 𝖥_1) + (∂_1(𝖥 - 𝖥_1))^2 = 0.\n',
    '    2 ∂_2(𝖥 - 𝖥_1) + ∂_1^2(𝖥 + 𝖥_1) + (∂_1(𝖥 - 𝖥_1))^2 = 0.\n\nLEMMA_A_PROOF.md s.4 derives (4.2) from Murnaghan-Nakayama and the Frobenius formula; no bilinear identity or vertex-operator theorem is cited.\n')

# ---- s.7.8 Step 3 (shrink: LEMMA_A Lemma 3.2)
span('*Step 3 (the dilaton point).* [P, given string and dilaton]',
     '\n\n    G_h(y) = w^{2h-2} Q_h(w - 1),',
     '*Step 3 (the dilaton point).* [P, given string and dilaton: LEMMA_A_PROOF.md Lemma 3.2] For h >= 2 put Q_h(s) := Sum_k (s^k/k!) ⟨Omega (ψ^2/(1-ψ))^k⟩. Lemma 3.2 there, eq. (3.1) at U = ∅, J = 0 (its G_U includes the a = 0 term), gives')
span('Genus 1: with a marked leg 1/(1-ψ_0),',
     'hence DG_1 = yH = u w Q_1(s; (1)),',
     'Genus 1: put H(y) := Sum_a M_1(1^{a+1}) y^a/a!. Lemma 3.2 there at U = (1), J = 1 gives H = w e^u Q_1(s; (1)), ')
rep('Every x-coefficient of D𝖥 is a Laurent polynomial in w, and Laurent polynomials embed in Q(τ)[[y]] because w = 1 + y + ... is non-constant. So (*) holds coefficientwise in the monomials x^a w^b.',
    'Every x-coefficient of D𝖥 and, modulo constants in y, of ∂_τ|_y 𝖥 is a Laurent polynomial in w: w does not depend on τ at fixed y, so ∂_τ|_y only differentiates the coefficients, and the log v of G_1 has a τ-independent coefficient. Laurent polynomials embed in Q(τ)[[y]], because w = 1 + y + ... is transcendental over Q(τ). So (*), in which constants are killed by D or by (1 - e^{βD}), holds coefficientwise in the monomials x^a w^b.')

# ---- s.7.8 Step 4 (D_VI)
rep('Each coefficient has a unique preimage, and N-powers are bounded above at fixed ε-order.',
    'Each coefficient has a unique preimage, and N-powers are bounded above at fixed ε-order (this is for D𝖥 and ∂_τ|_y 𝖥 with the constants in y removed; a constant x^{2h-2} G_h(0) would map to ε^{2-2h} N^0, which is why they must be removed first).')
rep('one has 𝓛[-N^2/(4σ)] = 0 identically, so the coefficient of N^2 is unconstrained by the next equation, and 𝓛[-(1/24) log N]',
    'one has 𝓛[c(τ) N^2] = 0 if and only if c = C/σ with C constant, so (E1) does not see the constant C = -1/4 (its value comes from Step 3, G_0, and only the fact that it has this form is used here); and 𝓛[-(1/24) log N]')

# ---- s.7.8 Status
rep('Status. One referee pass (logic, literature and independent computation in one lens: sound with minor fixes, all applied above). The referee confirmed the value by an MV-free computation (Mumford/Faber recursion and DVV) at h <= 4 at three τ and at h = 5 at one τ [E]. Scope: Step 1 settles (B1) for all g and μ, the same derivation as Proposition MV(A) of the MVmatch report; it does not by itself',
    'Status. Two referee passes (round 5: logic, literature and independent computation in one lens; round 6: D_VI; both sound with minor fixes, all applied above). Round 6: (*), the order ε^0, ε^1 expansion and (E1) were re-derived by hand, and compared with Schur data at N^{-k}, k <= 14, at τ = 2/7, -1/3, 3, and k <= 16 at τ = 5/3, evidence only. The value was confirmed by an MV-free computation at h <= 9 at six τ [E] (round 6, Mumford/Faber + DVV, no Schur or MV input). Scope: (B1) for all g and μ is LEMMA_A_PROOF.md Proposition 2.1, the same derivation as Proposition MV(A) of the MVmatch report; it does not by itself')

# ---- s.8.1
rep('    Prop 2.10 (symplecticity, K = 1, orthogonality)  P           Theorem A only; NEW in round 5, one referee pass\n',
    '    Prop 2.10 (symplecticity, K = 1, orthogonality)  P           Theorem A only; round 5, re-derived by referees in rounds 5 and 6\n')
rep('    Prop 2.11 (Givental-Teleman on T_F)              P mod C     (C1)-(C3), (C5); not used for Theorems LC_P2, R_P2\n',
    '    Prop 2.11 (Givental-Teleman on T_F)              P mod C     (C1)-(C5); not used for Theorems LC_P2, R_P2; not refereed in round 6\n')
rep('    Prop 2.12 (lem:modrel, both relations)           P           no citation; NEW in round 5, one referee pass;\n',
    '    Prop 2.12 (lem:modrel, both relations)           P           no citation; round 5, re-derived by referees in rounds 5 and 6;\n')
rep('    Lemma R, Theorem E00                             P           s.3.1, H2 Steps 1-2, Prop 2.10\n',
    '    Lemma R, Theorem E00                             P           s.3.1, H2 Steps 1-2 (and Step 4 for E00(iv)), Prop 2.10\n')
rep('    Fpot_h on 𝔇, (H0)                                P           dilaton equation; NEW in this edit, unrefereed\n',
    '    Fpot_h on 𝔇, (H0)                                P           dilaton equation; round 5, re-derived line by line in round 6 (D_rows, D_hodge), sound\n')
span('    Lemma A_P2(b) (vanishing)                        C ',
     '    Theorem VI (vertex identity (VI))',
     '    Lemma A_P2(b) (vanishing)                        P mod (MV)  lemmaA/LEMMA_A_PROOF.md Theorem L, and textbook facts; its open point O1\n'
     '                                                                 answered by the round-6 referee L_logic (no major item);\n'
     '                                                                 see s.7.4 and 8.2\n'
     '    Theorem R_P2                                     P mod (MV)  Theorem LC_P2, (H0), (S_P2), (C_P2), Lemma 7.1, Lemma A_P2(b) for h <= g-1\n')
span('    Theorem VI (vertex identity (VI))                P mod C ',
     '    Proposition 13.2_P2 (G_g = GW^cf_g)',
     '    Theorem VI (vertex identity (VI))                P mod (MV)  LLZ eq. (1) only, at μ = 1^a; textbook facts (string, dilaton, Frobenius,\n'
     '                                                                 Murnaghan-Nakayama, principal specialisation); (7.2);\n'
     '                                                                 two referee passes (rounds 5, 6)\n'
     '    Corollary SV_P2                                  P mod (MV)  above + (H0) + Lemma A_P2(b) at h = g, Theorem VI\n')

# ---- s.8.2
rep('The following is the outcome of a dependency reading, not of a full referee pass over the document.',
    'The following is the outcome of a dependency reading, confirmed for s.3-s.7 by the round-6 referees D_rows, D_core and D_hodge.')
rep('LLZ eq. (1) (arXiv v3 read, JDG numbering not confirmed), with OP Theorem 1',
    'LLZ eq. (1) (JDG 65 (2003) and arXiv v3 read, same numbering), with OP Theorem 1')
rep('Open there: O1, the assembled text has not had its own referee; O2, citation details, none mathematical; O3, the proof of (MV) in LLZ was not re-checked.',
    'Its O1 (own referee) is answered by the round-6 pass L_logic (sound with minor fixes, no major item). Open there: O2, citation details, none mathematical; O3, the proof of (MV) in LLZ was not re-checked.')
rep('LLZ eq. (1) and Theorem 1 eq. (5); the boson-fermion correspondence (vertex-operator form of ψ, ψ^*); the principal specialisation. The string and dilaton equations and (7.2) are standard or proved here.',
    'LLZ eq. (1), only at μ = 1^a, through LEMMA_A_PROOF.md Proposition 2.1; that file\'s Lemma 3.2 and eq. (4.2), whose proofs use neither (MV) nor Theorem L; textbook facts (string, dilaton, Frobenius, Murnaghan-Nakayama, principal specialisation). (7.2) is proved here.')
rep('Unconfirmed citations are collected in s.7.4 and s.7.8: JDG numbering of LLZ; theorem numbers for the mKP bilinear identity in Date-Jimbo-Kashiwara-Miwa and Miwa-Jimbo-Date; Macdonald I.3 Ex. 2. F2 has',
    'Unconfirmed citations: Macdonald I.3 Ex. 2 (s.7.8); no textbook number is pinned for Frobenius\' transposition identity. F2 has')

# ---- s.9 item 8 (D_rows)
rep('8. **sec:cfgen (5140).** Cf is described through τ -> -1/(3τ), which does not fix the Hauptmodul y. Cf must be defined as the identity on G_{T_F} (BS\'s commented-out identification Y^cf = Y) and extended to A^{-1}. Here this is a stated convention (s.1.3).',
    '8. **sec:cfgen (5140).** Cf is defined on generators, Cf(A) = A_cf, Cf(F) = F_cf, so it is the identity on G_{T_F}, but BS adds the gloss Cf(f)(τ) = f(-1/(3τ)), which cannot be literal on G_{T_F} since it does not fix the Hauptmodul y. Read Cf through the generators (compare BS\'s commented-out identification Y^cf = Y) and extend it to A^{-1} (s.1.3).')

# ---- Appendix (D_core)
rep('and gap_identity_check.pl, dexrel_check.pl and the data comparisons are demoted to evidence.\n',
    'and gap_identity_check.pl, dexrel_check.pl and the data comparisons are demoted to evidence.\n\n'
    '**Caveat (E-scaling, s.3).** By the scaling law of s.3, row and edge factors at different weight pairs differ only by ε1 ε2 and the scale E. So row and edge checks "at n weight pairs" below are one test of the scale-free object plus a test of the scaling; they do not probe non-special framing. Checks of Hodge-side quantities do depend on τ.\n\n'
    '**Round 6 (SP/lp2_round6/: D_rows, D_core, D_hodge, D_VI, L_logic).** Referees\' fresh scripts, listed with their outcomes in the round-6 referee reports. Evidence only.\n')

if fails:
    print("FAILED", len(fails))
    for c, o in fails:
        print(c, repr(o))
    sys.exit(1)
open(F, "w", encoding="utf-8", newline="").write(s)
print("ok")
