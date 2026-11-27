# One notation for the paper

Working file for the section writers; not part of the paper. Every recurring symbol has a macro in `main.tex`. Source abbreviations: FB = FRAMEB_PROOF.md, DR = LP2_REGULARITY_PROOF_draft.md, LA = LEMMA_A_PROOF.md, BS = BS_v2_2026-11-02/main.tex, CP = leaf_structure_proof.md (conifold).

Rule of precedence: where BS has a notation, it is used. Otherwise the DR or FB letter is kept whenever it does not clash, because transcription errors are the main risk. A symbol with a clash is renamed according to the tables below. A local symbol not listed here (a dummy index, an auxiliary series inside one proof) keeps its source letter, provided it does not hit a symbol of s.1-s.6 below; if it does, rename it at the top of the proof and say so.

## 1. Weights, tori, target

| Paper | Macro | Meaning | FB | DR | LA | BS |
|---|---|---|---|---|---|---|
| ε1, ε2 | `\epsilon_1, \epsilon_2` | fibre weights | eps1, eps2 | ε1, ε2 | | same |
| E | `E` | ε1 + ε2, with ε1 ε2 E ≠ 0 | E | E | E (Cor L3) | written ε1+ε2 |
| E(y) | `\Ey` | BS's function; always with its argument | E(y) | E(y) | | E(y) |
| 𝔼 | `\Hdg` | Hodge bundle | 𝔼 | E_Hodge | | |
| Z | `Z` | K_P2 x C^2, line class H | Z | | | Z |
| Tmax, TF | `\Tmax`, `\TF` | maximal torus; fibrewise torus α_i = E/3 | T, T_max, T_F | T_F | | same |
| α* | `\alpha^*` | (E/3)(1,1,1) | α* | α* | | |
| wt_i | `\wts_i` | {3α_i, α_j-α_i, α_k-α_i, -ε1, -ε2} | wt_i | | | w_i^{(j)} |
| w_i (sans) | `\wprod_i` | product of wt_i | W_i^0 | | | 𝗐_i (BS's symbol) |
| **w** | `\bw` | (E, -ε1, -ε2) | w | w | | |
| **w**_τ | `\bwt` | (1, τ, -1-τ) | | | (1,τ,-1-τ) | |
| τ | `\tau` | framing -ε1/E | | τ | τ | BS's modular τ is written `\tau_{\bbH}`, only in a remark |
| σ_τ | `\st` | τ(1+τ); `\st'` = 1+2τ | | σ_τ, σ | sigma | |
| 𝒩 | `\Ntw` | O(-3)_{3α} + O_{-ε1} + O_{-ε2} | N | | | |
| e_k(α) | `e_k(\alpha)` | elementary symmetric; never bare e_k | e_k | e_k(α) | | same |
| H, H' | `H`, `H'` | equivariant line class; non-equivariant hyperplane class on TF | same | | | H |
| P_i | `P_i` | fixed points | same | | | same |
| 𝒬 | `\cQ` | Novikov and mirror variable | Q | 𝒬 | | 𝒬 |

## 2. Base, mirror data, conifold

| Paper | Macro | Meaning | Sources |
|---|---|---|---|
| θ | `\theta` | θ_y = y d/dy (always y) | FB theta, DR θ, BS θ_y |
| θ_𝒬 | `\thQ` | 𝒬 d/d𝒬 | FB theta_Q |
| ξ | `\xi` | 1 + 27y | FB Delta_c = xi; DR Δ_c; BS ξ. **Δ_c does not appear in the paper.** |
| X | `X` | 1/ξ | same in all |
| A, D, F | `A, D, F` | BS eq:ADEFXdef; F = θ log A | same |
| 𝖥, 𝖠 | `\sfF`, `\sfA` | indeterminates for F, A | FB 𝖥, 𝖠 |
| a_d | `a_d` | (3d)!/(d!)^3 | same |
| H_d | `H_d` | harmonic numbers | same |
| L̂(Y) | `\Lhat(Y)` | ∏(Y+α_j) + 27yY^3 | FB Lhat, DR L^, BS \widehat\LL |
| W(Y), W_i | `W(Y), W_i` | e_1Y²+2e_2Y+3e_3; W_i = W(Y_i) | same; on TF the sign is + (s.8 below) |
| Y_i, u_i | `Y_i`, `\mathsf{u}_i` | root; canonical coordinate θ𝗎_i = Y_i | FB u_j; BS 𝗎_i. FB's U_j stays `U_j` |
| t, t_i | `t`, `t_i` | y^{1/3}, ω^i t | same |
| ω | `\omega` | e^{2πi/3} | same. The TFT is `\TFT`, see s.4 |
| q_j | `q_j` | Y_j + α_j, critical coordinates | DR q_i |
| v_i, V_i | `v_i`, `V_i` | 3/(E+3Y_i); E v_i | same |
| 𝓛_2 | `\Ltwo` | θ²+3y(3θ+1)(3θ+2) | L2b |
| 𝓛_PF | `\LPF` | operator of eq:pf | L_pf |
| 𝒫 | `\Pconj` | e^{-u/z} 𝓛_PF e^{u/z} | DR, FB 𝒫 |
| w_cf | `\wcf` | solution of 𝓛_2 analytic at ξ=0, value 1 | wc |
| t_c | `\tc` | θ t_c = -w_cf, t_c(ξ=0) = 0 | t_c |
| c_cf | `\ccf` | θ t_c | C1cf; also DR's g of Lemma 4.2 |
| c_1 | `\cone` | non-vanishing section; c_1^{-1} replaces C1i | c1, C1i |
| F_c | `\Fc` | θ log c_1 | XX |
| d_2 | `\dtwo` | E(c_1^{-2}-1)/(2ξ) | d2 |
| t_cf, A_cf, F_cf | `\tcf, \Acf, \Fcf` | BS sec:cfgen | same |
| 𝗀, 𝗀_2 | `\gI`, `\gtwo` | I_0^{[1]} - log y; series part of I_1^{[2]} | FB g, g2 |
| 𝔥 | `\hI` | Σ 3a_d(H_{3d}-H_d)(-y)^d | DR, FB h |
| υ | `\upsilon` | -27y in the Barnes computation | FB u |
| l | `l` | log(-y) | same |
| ρ | `\rhoc` or `\rho` | unit -E v_0/t_c | same. FB's reduction map ρ is `\redn` |
| J, S | `\Jsm`, `\Fsol(H^a)` | small J, fundamental solution | BS 𝖩, 𝖲(H^a); FB S(a) |
| **t**, **t**_a | `\ft`, `\ft_a` | point of the state space, flat coordinates | FB t, t_i; BS t |
| ∘ | `\qprod` | quantum product; H∘ for FB's H* | BS ∘_t |
| 𝔢_i, 𝔣_i | `\eid_i`, `\fid_i` | idempotents, 𝔣_i = Δ_i^{1/2} 𝔢_i | BS; FB, DR e_i, f_i. **FB's e_i ⇒ 𝔢_i.** Elementary symmetric functions stay e_k(α) |
| Δ_i | `\Delta_i` | η(𝔢_i,𝔢_i)^{-1} | same |
| η | `\eta` | pairing | same |
| Ψ | `\Psi` | matrix of the 𝔣_i | same |
| 𝒲, dΩ | `\cW`, `\dd\Omega` | superpotential, primitive form | BS eq:spotP2 |
| 𝗉_𝒲 | `\Pcr` | y^{-1/3} e^{Σx/3} | DR P |
| ζ | `\zeta_j` | jet coordinates at a critical point | FB, DR ξ. **ξ means 1+27y only.** |
| 𝖧 | `\Hess` | Hessian, FB's H = diag(α_j+Y) - (Y/3)J, DR's H^{ab}; FB's G is written -𝖧^{-1} | H, G |
| Jet_m | `\Jet_m(\zeta)` | Σ(α_j+Y)ζ_j^m/m! - 3Y(Σζ_j/3)^m/m! | FB T_m |
| Asym, 𝒜, 𝒜_i | `\Asym`, `\Asy`, `\Asy_i` | formal Laplace expansion; 𝒜_i = c_0 W_i^{-1/2} ι_i(Φ) | DR 𝒜; FB A_i |
| Θ, h_v, h_w | `\Theta`, `h_{\mathsf v}`, `\hb` | BS Θ; h_w = Υ(z/ε1, z/ε2) of BS | same |
| N(z) | `\Nz(z)` | sqrt(ε1ε2) Θ(-z/ε1) Θ(-z/ε2) | N |
| Λ^∨, Hg | `\Lamv(u)`, `\Hg_{\mathsf v}`, `\Hg_{\bw}` | Hg_v = ∏Λ^∨(v_a)/v_a^g | FB Hodge_v; DR Omega_w; LA Omega_h = Hg at **w**_τ. **Omega with a subscript w or h does not appear.** |
| Ω_{g,n} | `\Omega_{g,n}` | the CohFT of Z | FB Omega_{g,n} |
| ω (TFT) | `\TFT` | ω(𝔢_i^n) = Δ_i^{g-1} | FB ω |
| Ω_∘ | `\Omega_{\circ}` | rank-one CohFT of C Step 1 | FB Theta_0 |

## 3. Fields and rings (FB renames)

| Paper | Macro | Definition | FB |
|---|---|---|---|
| 𝕂 | `\bbK` | Q(ε1,ε2) | K |
| 𝕂_y | `\Ky` | 𝕂(y) | K_y |
| 𝕃 | `\bbL` | 𝕂(t), t³=y | L |
| 𝕃_ω | `\Lom` | 𝕂(ω)(t) | L_omega |
| 𝕂_α | `\Kal` | Q(α,ε) | Kw |
| 𝕂̃_α | `\Kalt` | 𝕂_α((𝗐_i)^{1/2}) | K~w |
| 𝕂̂_α | `\Kalh` | the larger field of the oscillating-integral argument | K^ |
| 𝔽_∘ | `\Fzero` | 𝕂̃_α((𝒬)) | F_0 |
| 𝔸 | `\ringA` | Q[α,ε,y] | 𝔸 |
| 𝔮 | `\primeq` | prime (α_i - E/3) | 𝔮 |
| 𝕆 | `\locO` | 𝔸_𝔮; 𝕆/𝔮𝕆 = 𝕂_y | O |
| 𝕆_α | `\locOal` | Q[α,ε] localised at 𝔮 ∩ Q[α,ε] | O_w |
| 𝒞 | `\ringC` | 𝕆[Y]/(L̂) | 𝒞 |
| red | `\redn` | 𝕆_α[[y]] → 𝕂[[y]], α_i → E/3 coefficientwise | rho |
| π_𝒞 | `\redC` | 𝒞 → 𝕃, Y → -(E/3)/(1+3t) | pi |
| ι_j | `\iota_j` | 𝒞 → 𝕂_α((y)), Y → Y_j | iota_j |
| 𝕃' | `\bbL'` | 𝕃_ω(κ_∘)(c_1,𝖥), κ_∘² = 3ε1ε2E | L', kappa_0 |
| 𝕃̃, **G** | `\Lt`, `\Gring` | as in Definition (Cf) | L~, 𝐆 |
| ev_A, ev_cf | `\evA`, `\evcf` | evaluation maps | same |
| R_ω | `\mathscr{O}_{\omega}` | 𝕂(ω)[y^{±1},ξ^{-1}] | R_omega |
| 𝐑 | `\RRing`, `\RRing_0` | 𝐑_0[[t_c]] | DR 𝐑 |
| k | `\Bbbk` | constant field of the symplecticity proof | DR 𝔽 |

## 4. Givental data

| Paper | Macro | Meaning | Sources |
|---|---|---|---|
| R, R_{aj} | `\Rm`, `\Rm_{aj}` | R-matrix; BS 𝖱_{aj}; R_m → `\Rm_m` | FB R |
| R̃_{aj} | `\Rt_{aj}` | stripped flat rows e^{-u_j/z} 𝗌_{aj}; equal to BS's 𝖲_{aj} in eq:canR | FB R~ |
| 𝗌_{aj} | `\mathsf{s}_{aj}` | (𝗐_j)^{-1/2} 𝖲(H^a)|_{P_j} | FB s_{aj}; hat form `\widehat{\mathsf{s}}_{aj}` |
| **B** | `\Rinv(\psi)` | R^{-1}(ψ) | FB B(psi) |
| R' | `\RB` | frame-B R-matrix; R_{TF} is `\Rm_{\TF}` | FB R' |
| T_R(z) | `T_{\Rm}(z)` | z(1 - R^{-1}(z)1); T_{R'} for FB's T'; T_{R'}^{(i)} the coefficient of 𝔢_i | FB T(z), T', T^{(i)}. Action by translation by a series b: `T_{b}` as in PPZ; R-action by scalar h: `\Rm_{h}` |
| 𝖦_Ψ | `\GPsi` | Ψ^{-1}θ_𝒬Ψ | FB W_Psi |
| **Λ** | `\mathbf{\Lambda}` | diag of eigenvalues of H∘ | FB Lambda |
| L_i | `\mathsf{L}_i` | eigenvalue of H∘ with value -α_i at 𝒬=0 | FB lambda_j. **λ means Hodge classes.** |
| **M**(ψ) | `\Mdef(\psi)` | defect of the flatness proof; 𝖣_i, 𝖣 for FB's D_i, D | FB M |
| 𝗌𝗁 | `\shQ` | ring map 𝒬 → 𝒬e^{**t**_1} | FB sigma |
| 𝗄, 𝗄_0 | `\mathsf{k}, \mathsf{k}_0` | constants with 𝗄²=-1 and 𝗄_0²=-1/3 | FB, DR c, c_0 |
| X_h | `\mathsf{X}_h` | Σ b_k((2k)!/B_{2k}) ch_{2k-1}(𝔼) | FB X_h |
| K_π | `K_{\pi}` | c_1(ω_π) in the marked Mumford formula | FB K |
| 𝐈_3 | `\mathrm{I}_3` | ∫λ_1λ_2λ_3 at g=3 | FB I |

## 5. Rows and factors

| Paper | Macro | Meaning | Sources |
|---|---|---|---|
| S(v,z), S^{[m]} | `\Srow(v,z)`, `\Srow^{[m]}` | closed form | same. Do not confuse with `\Fsol` |
| S_∞ | `\Sinf` | | same |
| r_l | `r_l` | recursion (Rec) | same |
| Φ, 𝖯_0 | `\Phi`, `\Pn_0` | algebraic rows, 𝖯_0 = Θ(z/ε1)Θ(z/ε2)Φ | FB Phi, P_0 |
| 𝖯_a^{(i)} | `\Pn_a^{(i)}` | Δ_i^{1/2} R_{ai} at generic weights | FB P_a^{(i)} |
| 𝖯̄_a^{(i)} | `\Pb_a^{(i)}` | TF rows, point frame | FB Pbar |
| 𝖯̂_a^{(i)} | `\Ph_a^{(i)}` | normalised frame-B rows | FB Phat, DR P^ |
| 𝖯^B_a | `\PB_a^{(i)}` | unnormalised frame-B rows | DR P_a^{B,(i)} |
| R^LR_a | `\RLR_a^{(i)}` | h_w 𝖯^B_a | DR P_a^{(i)} (large-radius frame). **DR's P_a^{(i)} is not FB's.** |
| 𝔻, 𝔻̄, 𝔻_i | `\opD`, `\opDb`, `\opD_i` | zθ + Y - (z/2)θ log W; the TF form with θ log Δ_i | FB 𝐓, 𝐓bar, T_i |
| Th_i | `\Thp_i` | Y_i + zθ | DR Th_i |
| **Q**, **Q̄**, **Q̂**_i | `\Qa`, `\Qb`, `\Qh_i` | 𝔻𝖯_0, 𝔻̄𝖯̄_0, 𝔻_i S(v_i,z) | FB 𝐐, Qbar, Q_i |
| 𝖰_1 | `\mathsf{Q}_1` | Th_i 𝖯^B_0 | DR Q_1 |
| **σ** | `\brow` | ρ𝖯̂^{(0)}(-ψ), components σ_a, coefficients σ^{(k)} | DR σ |
| Ŷ | `\Yh` | c_1^{-1} Y | DR Ŷ |
| Dil, Dil̄, Dil̂ | `\Dil_i(d)`, `\Dilb_i(d)`, `\Dilh_i(d)` | dilaton factors in the three frames | FB Dil, Dilbar, DR Dil^. BS's 𝖣_{i,n}: Δ_i^{1/2} 𝖣_{i,n} = Dil̄_i(n) |
| 𝓔, 𝓔̄, 𝓔̂ | `\Edge_{ij}`, `\Edgeb_{ij}`, `\Edgeh_{ij}` | edge factors in the three frames | FB E_ij, Ebar, DR Ehat. BS's 𝖤^{ij} times Δ_i^{1/2}Δ_j^{1/2} is 𝓔̄ |
| 𝓔̃ | `\Edget` | t_c²ρ²𝓔̂_00 | DR Ebar in Theorem E00. **Not FB's Ebar.** |
| 𝒫_{ij} | `\cP_{ij}` | -⟨𝖯̂^{(i)}(-ψ), 𝖯̂^{(j)}(-w)⟩ | DR P_ij |
| 𝐯_{j,k} | `\mathbf{v}_{j,k}` | [t_c^k] 𝖯̂^{(j)} | DR v_{j,k} |
| 𝖢(ψ,ψ') | `\mathsf{C}(\psi,\psi')` | the series of the residue lemma | DR P(ψ,ψ') |
| 𝒜_0 | `\cA_0` | matrix of the closure statement, if kept | DR A_0 |

## 6. Graphs and sums

| Paper | Macro | Meaning | Sources |
|---|---|---|---|
| Ver(Γ), Edg(Γ) | `\Ver(\Gamma)`, `\Edg(\Gamma)` | vertices, edges | FB, DR V, E; BS 𝒱_Γ, ℰ_Γ |
| 𝔠 | `\col`, `\col_v` | colouring Ver(Γ) → {0,1,2} | c, c_v. Colour 0 is the singular colour on the branch of t with t → -1/3 |
| g_v, n_v, m_v, a_h, d_k | same | as in (0.1) | same |
| Cont(Γ,𝔠) | `\GCont(\Gamma,\col)` | decoration sum of one coloured graph | DR T(Γ,c); BS Cont |
| 𝒢_g, 𝒢_g[𝖥] | `\Gsum_g`, `\Gsum_g[\sfF]` | Hodge-vertex frame-B sum | G_g, G_g[𝖥] |
| 𝒯_g, 𝒯_{Γ,dec} | `\Tsum_g`, `\Tsum_{\Gamma,\mathrm{dec}}` | polynomial over 𝕆 | FB 𝒯_g, T_{Γ,dec} |
| 𝒯̄_g[𝖥], 𝒯̄_g(F) | `\Tsumb_g[\sfF]` | point-frame sum on TF; 𝒯̄_g(F) is its value at 𝖥=F | FB Tbar_g[𝖥], 𝒯bar_g(F) |
| ℱ_g^{Tmax}, ℱ_g^{TF} | `\sF_g^{\Tmax}`, `\sF_g^{\TF}` | Σ_{d≥0} 𝒬^d GW_{g,d}, including d=0 | FB F^T_g, F_g^{T_F}. BS's symbol is GW_g with d>0; see s.8 |
| ℱ_{k,h}, ℱ^cf_{k,h} | `\sF_{k,h}`, `\sFcf_{k,h}` | coefficients; BS ℱ_{k,g} with g written h | same |
| GW^cf_g | `\GWcf_g` | Cf(ℱ_g^{TF}) | same |
| GW_{g,d} | `\GW_{g,d}^{\Tmax}` | invariants | FB GW^T_{g,d} |
| 𝔟_{k,h} | `\bfrak_{k,h}` | | same |
| β_n | `\beta_n(\epsilon_1,\epsilon_2)` | BS eq:bkgcf | same |

## 7. Leaves, Hodge side, Mariño-Vafa

| Paper | Macro | Meaning | Sources |
|---|---|---|---|
| **λ**_j | `\lf_j(\psi)` | 1/(E - jψ), **λ**_0 = 1/E | DR λ_j. **λ_i means Hodge classes only.** |
| **λ**^τ_j | `\lf^{\tau}_j` | 1/(1 - jψ), the scale-1 leaf | LA lam_j |
| 𝓛_m, 𝓛_m^+ | `\Lsp_m`, `\Lsp_m^{+}` | span of **λ**_1..**λ**_m; ℂ + 𝓛_m | DR Λ_m, Λ_m^+ |
| 𝗊_∘ | `\qo` | E - E²**λ**_1 | DR q_0. **q_0 means Y_0 + E/3.** |
| 𝗊_t, 𝗊_t^{(i)} | `\qt`, `\qt^{(i)}` | -ψ S(v_i,-ψ) | DR q_t |
| T_dil | `\Tdil` | 𝗊_∘ + ψ = -ψ²/(E-ψ); at scale 1 `\Tdil^{\tau}` = -ψ²/(1-ψ); `T^{+}` = -`\Tdil^{\tau}` | DR T, LA T, Tp |
| F^pot_h | `\Fpot_h` | | Fpot_h |
| 𝔡^s | `\dbase^s` | base-point derivatives | 𝔡^s |
| 𝔇 | `\Dom` | 𝗊_∘ + t_c 𝐑[[ψ]] | 𝔇 |
| 𝖬, 𝖴, 𝖵 | `\msM`, `\msU`, `\msV` | multisets: 𝖬 has integer parts ≥ 2 | DR, LA R; LA U, V (and T in old notes) |
| Ξ_h(𝖴) | `\Xiz_h(\msU)` | scale E | DR Xi_h |
| Ξ^τ_h(𝖴) | `\Xit_h(\msU)` | scale 1, class Hg at **w**_τ | LA Xi_h, X_g for g ≥ 1 |
| X_0(𝖴) | `\mathrm{X}_0(\msU)` | genus-0 formula, not a Ξ | LA X_0 |
| dg | `\dg(\msU)`, `\dg_v` | Σ(j-1) | D(U), D(R), D_v |
| δ_g(𝖴) | `\delta_g(\msU)` | dg(𝖴) - 2g + 2 | same |
| J_𝖴, r_𝖴 | `J_{\msU}`, `r_{\msU}` | | same |
| Z[h,n], Z_h(τ) | `\ZZ[h,n]`, `\ZZ_h(\tau)` | Ξ_h(1^n); value at scale 1 | same |
| M_irr(h) | `\Mirr(h)` | | M_irr(h) |
| x | `x` | MV expansion variable; LLZ's λ is `\lambda_{\mathrm{MV}}` with x = iλ_MV | same. Log coordinates are `x_j` with subscript |
| ℏ | `\hbar` | σ_τ x² | same |
| 𝗒, 𝗒_k | `\yv`, `\yv_k` | x σ_τ p_1; x σ_τ Π_k p_k | LA y, y_k. **y means the base coordinate.** |
| ϑ | `\thv` | 𝗒 d/d𝗒 | LA theta, DR-s.7.8 D. DR's formal variable ϑ in (H0) becomes `\mathsf{r}` |
| 𝗎 | `\tv` | 𝗒 = 𝗎 e^{-𝗎} | LA u |
| ϖ, ϖ̌ | `\pw`, `\pv` | 1/(1-𝗎); 1/ϖ | LA, DR w; DR v |
| x̌ | `\xinv` | 1/x | DR ε in s.7.8 |
| s, N, M, ℓ | `s, N, M, \ell` | LA s = ϖ - 1; DR s.7.8 N, M = N+τ, ℓ = log(M/N) | same |
| 𝖰_g(s;𝖴) | `\mathsf{Q}_g` | LA Q_g | same. No clash in s08-s09 with `\mathsf{Q}_1` of s03 because of the two arguments; write `\mathsf{Q}_g(s;\msU)` |
| 𝒵, ℱ_MV | `\MVZ`, `\MVF` | MV partition function, its log | LA Z, F; DR 𝖹, 𝖥. **𝖥 means the indeterminate only.** |
| 𝒵_0, 𝒵_1 | `\MVZ_0`, `\MVZ_1` | | LA Z_0, Z_1 |
| G_h, G_𝖴^{(g)}, M_g(μ) | `G_h(\yv)`, `G^{(g)}_{\msU}(\yv)`, `M_g(\mu)` | | same |
| 𝖠̃ etc. | `\mathfrak{A}, \mathfrak{B}, \mathfrak{C}` | grading pieces ℱ_MV = x𝔄 + 𝔅 + x̌ℭ + ... | DR 𝖠, 𝖡, 𝖢 |
| β_τ | `\gamma_{\tau}` | τx | LA beta |
| 𝖖 | `\mathsf{q}` | e^x, Schur specialisation `\mathsf{q}^{\bwy}` | LA q, ρ |
| ψ_i, ψ*_i | `\boldsymbol{\psi}_i`, `\boldsymbol{\psi}^*_i` | fermions in the derivation of (4.2) | LA psi_i. Digamma is `\psi^{(0)}`, `\psi^{(1)}` |
| a_k | `\mathrm{a}_k` | bosons | LA alpha_k |
| 𝗍_k | `\mathsf{t}_k` | p_k/k | LA t_k |
| ς_k | `\varsigma_k` | Hirota shift | LA eps_k |
| ℍ_k | `\Hir_k` | Hirota derivative | LA H_k |
| 𝖡_k | `\Bil_k(\Hir)` | bilinear polynomials | LA P_k |
| Sch_j, Eu | `\Sch_j`, `\Eul` | | same |
| 𝒢(ς) | `\cG(\varsigma)` | log-shift generating function | LA Phi(eps) |
| φ_B, B | `\varphi_B`, `B` | block factor, block | LA phi_B |
| ϱ | `\varrho` | partition of k_0 | LA alpha |
| 𝗋_i | `\mathsf{r}_i` | Chern roots of 𝔼 | LA alpha_i |
| 𝖪 | `\mathsf{K}` | content-weight operator | LA W |
| h(□) | `\mathsf{h}(\square)` | hook length | h |
| 𝗍 | `\mathsf{t}` | Borel variable in s.7.8 (no subscript) | DR t |
| f_Δ | `f_{\Delta}` | f_s - f_reg | DR W_0 |

## 8. Conventions where the paper's reading differs from BS (neutral wording; no errata list)

State each once, at the place named, as a convention or a brief remark. Sources: FB s.5.1, DR s.9.

- d ≥ 0: ℱ_g includes d = 0. BS's GW_g has d > 0. Stated in s02; constant in s05; effect on O(1) in s10.
- Double bracket ⟨⟨...⟩⟩ with Σ_{d≥0,m≥0} 𝒬^d/m!: s04, before the genus-0 differential equation.
- f_j^* is the constant covector at 𝒬 = 0: s04.
- W_i on TF is +(E/9)(3Y_i+E)²: s02. No mention of a misprint.
- Cf is the substitution on generators, identity on coefficients, with X and Y_i fixed: Definition in s06. Do not discuss the modular formula as an error. If a remark is wanted: "the substitution is not the pull-back under τ → -1/(3τ) on the whole ring".
- Branch log y = log(-y) + iπ for t_cf: Lemma in s06.
- Source of the mirror theorem at generic weights: [CCIT15, Theorem 31]; [Coa09] covers TF: s04, one sentence.
- Π_{γ_cf} read through its expansion: s10, statement of the G_{k,h} corollary.
- The constant 𝗄 with 𝗄² = -1 in the oscillating-integral normalisation: s04.
- Coefficients in Q[X^{±1},(1-X)^{-1}] rather than Q[X^{±1}]: s01 "not claimed" and s10.

## 9. BS references: phrase to use, label for \draftnote

Line numbers must not appear. Use the phrase, and if useful `\draftnote{BS label: ...}`.

| Source reference | Phrase | BS label |
|---|---|---|
| eq:Ifun, eq:pf, eq:Ian | the I-function of Z; the equivariant Picard-Fuchs equation; its expansion | eq:Ifun, eq:pf, eq:Ian |
| eq:ADEFXdef | the generators A, D, E(y), F, X | eq:ADEFXdef |
| eq:etaflat, eq:qprod, eq:fundsol, eq:Jred | the pairing; quantum product; fundamental solution; small J-function | same |
| eq:J=I, eq:SHb, eq:RHb, eq:canR | mirror theorem; components of the fundamental solution; row recursion; canonical R-matrix | same |
| line 1233, 1277 | the definition of the refined generating series by localisation | defn: refined GW series |
| line 2349, 2383 | weights of the maximal torus; tangent characters | no label |
| line 2795-2798, eq:tau | the modular parameter | eq:tau |
| eq:spotP2, line 4256, 4245-4252 | mirror superpotential; Asym | eq:spotP2 |
| exmpl:stirling, line 4278, 4282 | formal Stirling; Θ, Υ | as stated |
| eq:normNLR | the normalisation of the oscillating integral | eq:normNLR |
| eq:gtH, eq:leafterms, eq:edgeterm | the stable graph sum and its factors | same |
| eq:F2ref | the genus-2 generating series | eq:F2ref |
| sec:cfgen, line 5145, 5172 | generators at the conifold point | sec:cfgen |
| line 5340 | the formula for Y_0^cf in the proof of the conifold leading term | thm:cfasym |
| thm:cfasym, lem:cfpole, rmk:cfgap | the conifold leading term; pole order; the refined conifold gap conjecture | same |
| def:dirint, line 927-929, 937, 955-959 | the properties of refined direct integration; Conifold asymptotics | def:dirint, eq:bkgcf |
| eq:pfinf, line 2755-2760 | the leading term of the Picard-Fuchs equation at z = ∞ | eq:pfinf |
| prop:cfgraphsum, eq:betagiv | the graph sum for β | same |
| prop:algindABC, eq:XFABC, eq:ABCmod, eq:GTFS3 | the algebraic independence statement; the invariant ring | same |
