import io
p='patch_d.py'; s=open(p,encoding='utf8').read()
rep=[
("𝖲_1 f = -l f'(M), l := log(M/N)","𝖲_1 f = -ℓ f'(M), ℓ := log(M/N)"),
("- coefficient of l: 2𝖠_τ'(M)","- coefficient of ℓ: 2𝖠_τ'(M)"),
("(1/σ^2) { -l (MN - σ'/12)","(1/σ^2) { -ℓ (MN - σ'/12)"),
("Bor[(1 - e^{-τt})/t] = l. Set","Bor[(1 - e^{-τt})/t] = ℓ. Set"),
("MN l - τ(M + N)/2 = Bor","MN ℓ - τ(M + N)/2 = Bor"),
("Using ∂_τ(-1/σ) + 2(-1/σ)/a_2 ... (in the form σ'/σ^2 - 2τ/σ^2 = 1/σ^2) and ∂_τ d_0 = -σ'/(24σ^2), one gets","Using σ'/σ^2 - 2τ/σ^2 = 1/σ^2 for the coefficient of 𝖯/t^3, and ∂_τ d_0 = -σ'/(24σ^2), one gets"),
("Setting t = a + b, t' = a - b and b = 0 at the z^{-2} coefficient gives","Setting t = a + b, t' = a - b and taking the coefficient at b = 0 (the z^{-2} coefficient is (1/2)(D_1^2 - D_2) τ_1 · τ_0, D_k the Hirota derivatives) gives"),
("Round 1's scope wording for this result: Step 1 settles","Scope: Step 1 settles"),
("(7.2). The tau-symmetry of Omega fixes the integration constants. Labels: [P] proved here, [C] cited. The mKP equation is the form allowed by the rulings.","(7.2). The proof combines the MV formula, the first mKP equation (KP/Toda is allowed by the rulings) and the cut-and-join eigenvalue on the p_1-slice, extracts the lowest order at the dilaton point, and uses the τ-symmetry of Omega to fix the integration constants."),
("and this vanishes at two complex values of σ_τ. No proof step uses the non-vanishing.","and this vanishes at two real values σ_τ < -1/4, i.e. at complex τ. No proof step uses the non-vanishing."),
("""    𝔡^{s+1}(g_1, ..., g_s, q_0) = Sum_m (1/m!) [m - (2h - 2 + s + m)] 𝔡-term_m = (2 - 2h - s) 𝔡^s(g_1, ..., g_s),

where 𝔡-term_m is the summand with m insertions of T. The m-dependence cancels.""","""    𝔡^{s+1}(g_1, ..., g_s, q_0) = Sum_m [m - (2h - 2 + s + m)] A_m = (2 - 2h - s) 𝔡^s(g_1, ..., g_s),

where A_m := (1/m!) Int Omega_w prod g_i prod_{k<=m} T is the summand of 𝔡^s with m insertions of T. The m-dependence cancels."""),
]
for a,b in rep:
    assert s.count(a)==1,a
    s=s.replace(a,b)
open(p,'w',encoding='utf8').write(s)
p='patch_b.py'; s=open(p,encoding='utf8').read()
rep=[("It has had one referee pass (sound; the referee re-derived Theorem A and S0-S6 by hand).","It has had one referee pass (sound with minor fixes, applied here; the referee re-derived Theorem A and S0-S6 by hand)."),
("- (c) Where lem:modrel does enter BS downstream:","- (c) (In this remark e_k := e_k(α).) Where lem:modrel does enter BS downstream:"),
]
for a,b in rep:
    assert s.count(a)==1,a
    s=s.replace(a,b)
open(p,'w',encoding='utf8').write(s)
