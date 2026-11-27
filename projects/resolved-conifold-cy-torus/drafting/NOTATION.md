# One notation for the note

Working file for the section writers; not part of the note. Every recurring symbol has a macro in `main.tex`; a symbol used in one section only may be defined at the top of that section file (as `\newcommand` or `\providecommand`, prefixed with the section's label prefix if it could clash). Source abbreviations: LCY = cy_torus/LC_CY_PROOF.md, RIG = cy_torus/RIGIDITY_PROOF.md, LSP = leaf_proof/leaf_structure_proof.md, LAP = lemmaA/LEMMA_A_PROOF.md, LP2 = lp2_explore/LP2_REGULARITY_PROOF_draft.md, CP = the companion LaTeX paper (lp2_explore/paper/), BS = BS_v2_2026-11-02/main.tex.

Rules of precedence. (1) Where BS has a notation for the same object, it is used (weights `\alpha_i`, `\epsilon_i`; `\cQ`; `\TCY`; `\Lamv`; `\beta_n`; `\GW`). (2) Where the companion paper CP has a macro for the same object, the macro and its meaning are kept (leaves, Hodge class, graph sums, rows, factors, fields). (3) Otherwise the LCY/RIG letter is kept when it does not clash. (4) A clash is resolved by the tables below; a renamed symbol is never used under its old name. Dummy indices and auxiliary series inside one proof keep their source letters, unless they hit a symbol of sections 1-8 below; then rename at the top of the proof and say so.

Two frames exist in the sources (point frame, frame B). In the note only the frame-B factors carry names (`\Dil`, `\Edge`, `\Delta`); point-frame rows carry the superscript pt (`\rpt`) and the rows with the Q = 0 constant carry Z (`\rZ`). "The graph sum" means the frame-B Hodge-vertex sum unless "point frame" is said.

## 1. Weights, tori, target

| Note | Macro | Meaning | LCY/RIG | LSP | BS |
|---|---|---|---|---|---|
| alpha_0, alpha_1, alpha_2 | `\alpha_0,\alpha_1,\alpha_2` | base weights; tangent weights of the zero section and of O(-1)^2 at p_+ | a0, a1, a2 | alpha0, alpha1, alpha2 (s.13) ; x, 1 on T_F | alpha_i |
| epsilon_1, epsilon_2 | `\epsilon_1,\epsilon_2` | weights of the two trivial directions (C^2) | e1, e2 | eps1, eps2 | epsilon_1, epsilon_2 |
| E | `E` | alpha_0 + alpha_1 + alpha_2, equal to epsilon_1 + epsilon_2 on T_CY. On the five-torus E means alpha_0 + alpha_1 + alpha_2 (RIG's sig); say so where both are in play (s05, s07) | E | c (on T_F); sig (s.13) | epsilon_1 + epsilon_2 written out |
| A (sans) | `\wA` | alpha_0 + alpha_2 = E - alpha_1 | A | 1 (on T_F) | |
| B (sans) | `\wB` | alpha_0 + alpha_1 = E - alpha_2 | B | x (on T_F) | |
| u (sans) | `\uq` | (E - 2alpha_1)(E - 2alpha_2) = alpha_0^2 - (alpha_1 - alpha_2)^2; also 2alpha_0 E + 4alpha_1 alpha_2 = u + E^2 | u | -(1-x)^2 | |
| (ND) | text | alpha_1 alpha_2 A B epsilon_1 epsilon_2 E != 0 | (ND) | x != 0, -1 | |
| T | `\Tmax` | the five-torus (C^*)^5 | T | T | T (maximal torus) |
| T_CY | `\TCY` | the Calabi-Yau four-torus E = epsilon_1 + epsilon_2 | T_CY | | T_CY |
| T_3 | `\Tthree` | alpha_0 = 0 inside T_CY | T_3 | | |
| T_F | `\TF` | alpha_0 = 0, alpha_i = epsilon_i (fibrewise) | T_F | T_F | |
| X, Z | `X`, `Z` | resolved conifold; Z = X x C^2 | Z | Z | X, Z |
| p_+, p_- | `p_\pm` | fixed points of the zero section | p_+, p_- | p_+, p_- | |
| wt_+, wt_- | `\wts_\pm` | (alpha_0, alpha_1, alpha_2, -epsilon_1, -epsilon_2), (-alpha_0, B, A, -epsilon_1, -epsilon_2) | wt_+, wt_- | wt_+, wt_- | |
| e(wt) | `\eul(\wts_\pm)` | product of the five weights; e(wt_+) = alpha_0 alpha_1 alpha_2 epsilon_1 epsilon_2, e(wt_-) = -alpha_0 A B epsilon_1 epsilon_2 | e(wt) | e(wt) | |
| N_cf | `\Nb` | O(-1)_{alpha_1} + O(-1)_{alpha_2} + O_{-epsilon_1} + O_{-epsilon_2} | Nb | Nb | |
| **w** | `\bw` | (E, -epsilon_1, -epsilon_2), the frame-B class | w | w = (c, -x, -1) | |
| nw | `\nw` | (alpha_1, A, -epsilon_1, -epsilon_2) | nw | nw | |
| h_v | `h_{\mathsf{v}}`, `h_{(\ldots)}` | Exp[-Sum_k bk_k z^{2k-1} Sum_a v_a^{1-2k}]; h_{u+v} = h_u h_v, h_{-v} = h_v^{-1}, h_v(z) h_v(-z) = 1 | h_v | h_v | Theta(z/v) = h_{(v)}^{-1} |
| b_k (sans) | `\bk_k` | B_{2k}/(2k(2k-1)) | bk_k | bk_k | |
| h_w, h_nw, h_reg | `\hw`, `\hnw`, `\hreg` | h at **w**; at nw; h_{(alpha_1, alpha_2, A, B, -E, -E)} | h_w, h_nw, hreg | h_w, h_nw, h_w^{-2} | |
| Q | `\cQ` | Novikov variable, degree variable of P^1 | Q | Q | Q (calligraphic) |
| t | `t` | log Q; the conifold point is t = 0 | t | t | |
| tau | `\tau` | 1 - 1/Q = 1 - e^{-t}. **tau means this throughout.** The scale-1 parameter of [LP2] (-epsilon_1/E) is `\tauw`, and appears only in s06 when a result of [LP2] is quoted | tau | tau | |
| Disc(Q) | `\Disc(\cQ)` | alpha_0^2 + (2alpha_0(alpha_1 + alpha_2) + 4alpha_1 alpha_2) Q + (alpha_1 - alpha_2)^2 Q^2 | Disc | | |
| GW_{g,d} | `\GW_{g,d}` | genus-g, degree d[P^1] invariant on T_CY by localisation; superscript for another torus: `\GW^{\Tmax}_{g,d}` | GW_{g,d} | GW_{g,d} | GW_g(Z, G; beta) |
| F_g | `\sF_g` | Sum_{d>=0} GW_{g,d} Q^d, d = 0 included | F_g | F_g | GW_g(Z, G), d > 0 only |
| G_g (script) | `\sG_g` | F_g - F_g(0) | G | | |
| P_g | `\Pg`, `\Ph{h}` | [t^{2g}] epsilon_1 epsilon_2 t^2/(4 sinh(epsilon_1 t/2) sinh(epsilon_2 t/2)) | P_g | P_g | |
| beta_n | `\beta_n(\epsilon_1,\epsilon_2)` | (n-1)! [t^n] 1/(4 sinh(epsilon_1 t/2) sinh(epsilon_2 t/2)); beta_{2g-2} = (2g-3)! P_g/(epsilon_1 epsilon_2) | beta_n | beta_n | beta_n |
| gap_g | `\gapg` | beta_{2g-2}(epsilon_1, epsilon_2) | gap_g | gap_g | |
| S_o(y) | `\Ssh(y)` | y/(2 sinh(y/2)) | S(y) | S | |
| Li_{-n} | `\Li_{-n}(\cQ)` | Sum_d d^n Q^d | Li_{-n} | | Li |
| A_n (cal) | `\Aeul_n(\cQ)` | numerator: Li_{-n} = A_n/(1-Q)^{n+1} | A_n | | |

Dropped: RIG's K_w for Q(alpha_0, alpha_1, alpha_2, epsilon_1) is `\KCY`; RIG's K for Q(alpha_1, alpha_2, epsilon_1) is `\Kred`.

## 2. Fields and rings

| Note | Macro | Definition | Source |
|---|---|---|---|
| K_alpha | `\Kal` | Q(alpha_0, alpha_1, alpha_2, epsilon_1, epsilon_2), the five weights indeterminate | LCY Kw; LSP Kw |
| K_CY | `\KCY` | Q(alpha_0, alpha_1, alpha_2, epsilon_1), epsilon_2 := E - epsilon_1; function field of Lie(T_CY) | LCY K; RIG Kw |
| K (blackboard) | `\Kred` | Q(alpha_1, alpha_2, epsilon_1); K_CY = K(alpha_0) | RIG K |
| K_3 | `\Kthree` | Q(alpha_1, alpha_2, epsilon_1) as the function field of T_3 (same field as K, different role) | RIG K_3 |
| C((z)) | `\Cz` | formal Laurent series in z | LCY Kz |
| R_z, R_z^o | `\Rz`, `\Rzo` | C((z))((tau)); C[[z]]((tau)) | LCY Rz, Rz0 |
| O, M (script) | `\Oan`, `\Man` | convergent power series in t (= in tau); meromorphic germs at t = 0; O(D) analytic on a disc | LCY O, M |
| QA | `\QA` | K(Q)[U]/(P(U)), the quadratic algebra; its involution U -> U_1 + U_2 - U exchanges the idempotents | LCY QA, cconj |
| R (bold) | `\RRing` | K_CY[[t]] (or C[[t]]), coefficient ring of the Hodge side | RIG Rr |
| C (cal) | `\Ccov` | the double cover y^2 = Disc(Q) of P^1_Q; K_CY(C) = K_CY(Q)[U]/(P(U)) | RIG C |
| y (sans) | `\ydc` | y^2 = Disc(Q); 2BQU - (alpha_0 + Q(alpha_1 - alpha_2)) = +-y | RIG y |
| Q_1 | `\cQ_1` | the small root of Disc in K[[alpha_0]], of valuation 2 | RIG Q_1 |
| iota, refl, TRANSFER, (P), (T), (PF) | text | as in LSP s.1.4, with K_x -> K_alpha | LSP s.1.4 |

## 3. Genus-zero data

| Note | Macro | Meaning | LCY | LSP (T_F) |
|---|---|---|---|---|
| W(U) | `\cW(U)` | alpha_1 log(1 - QU) - alpha_2 log U + A log(1 - U); W_k := theta_U^k W | W, W_k | W, W_k |
| theta_U | `\thU` | U d/dU | theta | theta |
| theta_Q | `\thQ` | Q d/dQ = d/dt; a prime on functions of t is d/dt | ' | ' |
| P(U) (sans) | `\Pq(U)` | B Q U^2 - (alpha_0 + Q(alpha_1 - alpha_2)) U - alpha_2 (critical equation) | LeqU; RIG P(U) | LeqU |
| U_1, U_2 | `U_1`, `U_2` | roots; U_1 -> 1, U_2 -> -alpha_2/B as t -> 0; idempotent 1 is singular, 2 regular. At Q = 0 (s05) the roots are labelled U_+ -> -alpha_2/alpha_0, QU_- -> alpha_0/B; say which labelling is in force | U_1, U_2; U_+, U_- | same |
| N(s) (sans) | `\Ns(s)` | B tau s^2 - (E + alpha_0 tau) s + A, the critical polynomial in the Euler variable | N(s) | x tau s^2 - c s + 1 |
| s^*, s~ | `\sast`, `\stil` | roots of N; U_1 = 1 - tau s^*, U_2 = 1 - tau s~; s^*(0) = A/E | s*, s~ | s*, s~ |
| delta | `\delta` | the root of delta^2 = (E + alpha_0 tau)^2 - 4AB tau = E^2 - (E^2 + u) tau + alpha_0^2 tau^2 with delta(0) = E; delta = E + alpha_0 tau - 2B tau s^* | delta | delta |
| D, D', D'' | `\Ddisc`, `\Ddisc'`, `\Ddisc''` | disc of analyticity of delta; intersected with |tau| < 1; minus the origin | D, D', D'' | same |
| L(U), L_i | `\Lv(U)`, `\Lv_i` | -alpha_1 Q U/(1 - QU); eigenvalues of H o; L_1 = A/(tau s^*) - alpha_0 = E/t + l(t), L_2 = B s^* - alpha_0 | L(U), L_i | L_i = (xU_i + 1)/tau |
| S, Pi (sans) | `\Lsum`, `\Lprod` | L_1 + L_2 = E/tau - alpha_0; L_1 L_2 = alpha_1 alpha_2/tau (only if a name is needed) | Ssum, Pprod | |
| l(t) | `\lrem` | L_1 - E/t, in O; l(0) = E/2 - alpha_0 - alpha_1 alpha_2/E | l, l_0 | l |
| eta, <a,b> | `\eta`, `\langle a,b\rangle` | pairing in the flat basis {1, H}; eta^{-1} = epsilon_1 epsilon_2 [[alpha_0 alpha_1 alpha_2, alpha_1 alpha_2],[alpha_1 alpha_2, -E]]; covectors paired by eta^{-1} | eta, < , > | same |
| o | `\qprod` | quantum product; H o H = (E/tau - alpha_0) H - (alpha_1 alpha_2/tau) 1 | * | * |
| e_i, f_i (frak) | `\eid_i`, `\fid_i` | idempotents (H - L_j)/(L_i - L_j); normalised idempotents Delta_i^{1/2} e_i | e_i | e_i |
| Delta_i | `\Delta_i` | 1/eta(e_i, e_i) = epsilon_1 epsilon_2 (L_i - L_j)(alpha_1 alpha_2 - E L_i) = alpha_1 A epsilon_1 epsilon_2 W_2(U_i) | Delta_i | Delta_i |
| A (blackboard) | `\Aop` | [[0, 1],[-alpha_1 alpha_2/tau, E/tau - alpha_0]] on covectors (transpose of the matrix of H o). **Not** the weight A | A | A |
| chi_i (bold) | `\chiv_i` | (1, L_i), eigen-covectors of A | chi_i | chi_i |
| pi(a) (varpi) | `\pim(a)` | <a, chi_2> = epsilon_1 epsilon_2 B s^*(alpha_1 alpha_2 a_0 - tau L_2 a_1) | pi | pi |
| Theta | `\Thf` | Delta_2/(L_2 - L_1) = epsilon_1 epsilon_2 (alpha_1 alpha_2 - E L_2) = -epsilon_1 epsilon_2 B tau s^* L_2 | Theta | Theta = -x^3 tau s*^2 |
| mu, omega | `\muq`, `\omq` | E/tau - A; (E + z)/tau - (A + z) = mu + z(1 - tau)/tau | mu, omega | mu, omega |

## 4. Euler integral, Laplace calculus, Gauss (s03)

| Note | Macro | Meaning | LCY | LSP |
|---|---|---|---|---|
| s | `s` | Euler variable, U = 1 - tau s | s | s |
| phi_tau(s), phi_0 | `\phi_\tau(s)`, `\phi_0` | A log s + alpha_1 log(1 - s) - alpha_2 log(1 - tau s); phi_0''(A/E) = -E^3/(alpha_1 A) | phi_tau | phi_tau |
| phi-breve | `\dphi(\tau)` | phi_tau(s^*) - phi_0(A/E), in tau O; dphi' = A - B s^* = alpha_2 - L_2 | dphi | dphi |
| a (frak) | `\alg(\tau)` | (1 - tau s^*)(phi_tau''(s^*)/phi_0''(A/E))^{1/2}, the branch with value 1; alg^2 = A delta (E - B tau s^*)/(E^3 s^*) | alg | alg |
| rho | `\rho` | tau/(t alg) = (Delta_1 t^2/(-epsilon_1 epsilon_2 E^3))^{-1/2}, rho(0) = 1 | rho | rho |
| w_k (sans) | `\wrk_k` | tau^k W_k(U_1) (`\wr` is LaTeX's wreath product; do not use it) | wr_k | wr_k |
| Pol_n | `\Pol_n(p)` | Pol_0 = p, Pol_n = (p^2 - p) Pol_{n-1}' | Pol_n | Pol_n |
| A (script) | `\Acal` | ring of formal sums Sum K_{pq} v^p z^q with the weight p/2 + q | Acal | Acal |
| <<K>>_a | `\gauss{K}{a}` | Gaussian functional | <<K>>_a | same |
| N (script), N-hat, N^full | `\Lnorm[f,a;u_0]`, `\Lhat`, `\Lfull` | Gaussian-normalised formal Laplace expansion; normalised to leading coefficient 1; with the prefactor | Lnorm, Lhat, Lfull | same |
| Rr_f | `\Rr_f(v)` | Taylor remainder of order >= 3 of the phase | Rr_f | Rr_f |
| wick_m | `\wick_m` | the Wick polynomials in W_2^{-1}, W_3, ... | wickZ_m | wickZ_m |
| N (cal) | `\Nlap(z,\tau)` | Lhat[phi_tau, (1 - tau s)^{-1}; s^*] = U_1 Lnorm[...]; equals r^pt_{0,1} | N(z, tau) | N |
| A_beta | `\Abeta(z)` | Lnorm[phi_0, 1; A/E] = h_{(E, -alpha_1, -A)} E/(E + z) | Abeta | Abeta = h_w c/(c + z) |
| a, b, c (sans) | `\ga`, `\gb`, `\gc` | 1 + alpha_2/z, 1 + A/z, 2 + E/z, Gauss's parameters. **Not** alpha, beta, gamma (weights, beta_n) | alpha, beta, gamma | same |
| F (sans) | `\Fhyp(z,\tau)`, `\Fhyp_n` | 2F1(a, b; c; tau) = Sum F_n tau^n; F_n = prod_{j<=n} (alpha_2 + jz)(A + jz)/(n! z^n prod_{j=2}^{n+1} (E + jz)) | F, F_n | F, F_n |
| Gauss | `\Gauss[Y]` | tau(1 - tau) Y_tautau + (c - (a + b + 1) tau) Y_tau - ab Y | Gauss | Gauss |
| X (sans) | `\Ephi`, `\Ephi_n` | Exp[-dphi/z] in C[1/z][[tau]]; its tau^n coefficient | Ephi, X_n | Ephi, X_n |
| H_n, K_n (sans) | `\Hcoef_n`, `\Kcoef_n` | [tau^n](F X); [tau^n](F X alg) | H_n, Psi_n | H_n, Psi_n |
| k_E | `k_E` | s(1 - s)(1 - tau s)^{-2} in the amplitude identity | kE | kE |
| Dam | `\Dam` (= `\mathsf{D}_{\mathrm{am}}`, main.tex) | d_tau q + q (d_tau f)/z | Dam | Dam |
| m (sans) | `\amp` (= `\mathsf{m}`, main.tex) | amplitude of a Laplace expansion (s03, s04) | am | am |
| S_phi, Phi | `\Sph`, `\Phia` | eikonal phase; transport amplitude, F = Exp[S_phi/z] Phi | Sph, Phi | same |
| Eik, P_tr, R_tr, T_tr | `\Eik`, `\Ptr`, `\Rtr`, `\Ttr` | eikonal polynomial and transport coefficients; P_tr = delta, R_tr = delta (log alg)_tau | same | same |
| D_z, D (symbol) | `\Dz`, `\Dd` | z d/dt; abstract derivation with D a = a D + z a' | Dz, Dd | Dz, Dd |
| Pf | `\mathsf{Pf}` | (1 - tau)^{-alpha_1/z} tau^{E/z + 1} (conjugating factor, remark only) | Pf | Pf |
| kappa_ph | `\kappa_{\mathrm{ph}}` | -e^{alpha_1 log(-1)/z} (constant of the phase, used nowhere) | kappa_ph | kappa_ph |

## 5. Rows and factors (frame B)

| Note | Macro | Meaning | LCY | LSP |
|---|---|---|---|---|
| r^pt_{0,i} | `\rpt_{0,i}(z;t)` | Sum_m z^m wick_m(W_2(U_i), W_3(U_i), ...) | r^pt_{0,i} | r^pt_{0,i} |
| r^Z_{0,i} | `\rZ_{0,i}` | h_nw r^pt_{0,i} (the row with the Q = 0 constant; only in s05, s07) | h_nw r^pt (R_Z) | rows of rmatrixAtClosed |
| r_{a,i} | `\rr_{a,i}` | frame-B rows: r_{0,i} = h_{(alpha_1, A, -E)} r^pt_{0,i} = h_w^{-1} h_nw r^pt_{0,i}; r_{1,i} = (L_i + z D_i) r_{0,i} | r'_{a,i} = r_{a,i} | r_{a,i} = h_w^{-1} r^pt |
| r_i (bold) | `\vv_i(z)` | (r_{0,i}, r_{1,i}), a covector; r_i[p] := [z^p] r_i | v_i | v_i |
| D_i (sans) | `\Dop_i` | d/dt - Delta_i'/(2 Delta_i) | D_i | D_i |
| R, R_Z, R' | `\Rm`, `\RZ`, `\RB` | R-matrix at generic weights (Teleman); the matrix with rows r^Z; frame-B matrix R_Z h_w^{-1} | R, R_Z, R' | same |
| omega (bold) | `\TFT` | semisimple TFT omega(e_i^{(x)n}) = Delta_i^{g-1} | omega | omega |
| B (bold), M (bold) | `\Rinv(\psi)`, `\Mdef` | R^{-1}(psi); defect matrix (s05 only) | B, M | B, M |
| Dil_i(d) | `\Dil_i(d)` | (-1)^d r_{0,i}[d-1], d >= 1 (frame B) | Dil_i | Dil_i |
| q_t, q_t^{(i)} | `\qt`, `\qt^{(i)}` | Sum_d Dil_1(d) psi^d = -psi r_{0,1}(-psi; t); at idempotent 2 with superscript | q_t | q_t |
| P_ij (cal) | `\Pij_{ij}(z,w)` | -<r_i(-z), r_j(-w)> | P_ij | P_ij |
| E_ij (cal) | `\Edge_{ij}(z,w)`, `\Edge_{ij}(m,n)` | [P_ij(z,w) - P_ij(-w,w)]/(z + w); its coefficients. In leaf language E_11(psi, psi') means z = psi, w = psi' | E_ij, EdN | E_ij |
| f_{j,b} | `f_{j,b}(\psi)` | [w^b] E_1j(psi, w), the legs | f_{j,b} | f_{j,b} |
| K_1, K_2 | `\Kone`, `\Ktwo` | <r_i(z), r_i(-z)>/Delta_i | K_1, K_2 | K_1, K_2 |
| sigma (bold) | `\brow(\psi)` | rho r_1(-psi) = (sigma_0, sigma_1); sigma^{(k)} := [t^k] sigma | sigma | sigma |
| Y (bold) | `\Yv(\psi)` | t sigma(psi) | Y | Y |
| G, G_k (sans) | `\Gk(z;t)`, `\Gk_k` | rho r_{0,1} = (E/(E + z))(tau/t) F X (KEY FORMULA); [t^k] G | G, G_k | G, G_k |
| g, g_k, Dg | `\gl`, `\gl_k`, `\Dg` | G(-psi; t); [t^k] g; d/dt (t g) = Sum (k+1) g_k t^k | g, g_k, Dg | same |
| V_k (sans) | `\Rk_k(z)` | [t^k](z G) = z G_k | R_k | R_k |
| E-bar, N-bar | `\Ebar`, `\Nbar` | t^2 rho^2 E_11 = N-bar/(psi + psi'); N-bar = -epsilon_1 epsilon_2 E^3 - <Y(psi), Y(psi')> | Ebar, Nbar | Ebar, Nbar |
| a_j^{(k)}, g_j^{(k)} | `a_j^{(k)}`, `g_j^{(k)}` | leaf coefficients of [t^k](rho q_t) and of g_k; g_j^{(k)} = -(j/E) a_j^{(k)}; a_1^{(0)} = -E^2 | same | same |
| A_j(t), A-hat_j | `\Aleaf_j(t)`, `\Ahat_j` | Sum_k a_j^{(k)} t^k; t A_j | A_j, Ahat_j | same |
| e^{(k)}_{il}, e_{il}(t) | `e^{(k)}_{il}`, `e_{il}(t)` | loop coefficients of [t^k](rho^2 E_11) in lam_i(psi) lam_l(psi') | same | same |
| S_j^{(k)}, S_j(t) | `S_j^{(k)}`, `S_j(t)` | leaf coefficients of sigma^{(k)} (pairs) | same | same |
| p_k (frak), C_j (sans) | `\Psh_k`, `\Cj_j` | prod_{i=1}^k ((k+1) alpha_2 - iE)((k+1) A - iE); p_{j-1}/(j! (j-1)! E^{2j-4}) | P_k, C_j | P_k, C_j |
| zeta | `\zq(t)` | tau e^{dphi/E} = t + O(t^2) | zeta | zeta |
| d_{k,n} | `d_{k,n}` | [t^k]((tau/t) tau^n), rational | d_{k,n} | d_{k,n} |
| chi_12(w) | `\chi_{12}(w)` | -epsilon_1 epsilon_2 alpha_1 alpha_2 A B h_reg(w)^{-1}/(E + w), the t = 0 leg of E_12 | chi12 | chi12 |
| kappa_2 | `\kappa_2` | epsilon_1 epsilon_2 alpha_1 alpha_2 A B (table of values at t = 0) | kappa_2 | x^3 |
| f_2, f_0 | `f_2(\tau)`, `f_0(\tau)` | phi_tau''(s~), phi_tau(s~) at the regular root (s04, regular idempotent) | f2, f0 | f2, f0 |

## 6. Graphs and the graph sum

| Note | Macro | Meaning | Source |
|---|---|---|---|
| Gamma, Ver, Edg | `\Gamma`, `\Ver(\Gamma)`, `\Edg(\Gamma)` | connected stable graph of genus g without legs | RIG (3.1) |
| g_v, n_v, m_v, a_h, d_k | same | genus and valence of v; number of dilaton leaves at v; half-edge exponents; dilaton degrees | LSP s.0.1, (C10) |
| c (frak) | `\col`, `\col_v` | colouring Ver(Gamma) -> {1, 2}; colour 1 singular | RIG c |
| Cont(Gamma, c) | `\GCont(\Gamma,\col)` | decoration sum of one coloured graph: Sum_dec prod_v Delta_{c_v}^{g_v-1} Int Hg_w psi^a prod Dil prod E | RIG (3.1) |
| G_g (script) | `\Gsum_g` | Sum_Gamma |Aut Gamma|^{-1} Sum_c Cont(Gamma, c), the frame-B Hodge-vertex graph sum | "the graph sum" |
| Hg_w | `\Hg_{\bw}` | prod_{a=1}^3 Lamv_g(w_a)/w_a^g, the vertex class; Hg_v for any multiset | Omega_w |
| Lamv(u) | `\Lamv(u)` | Sum_i (-1)^i lambda_i u^{g-i} | Lambda^vee |
| E (blackboard) | `\Hdg` | Hodge bundle | E |

## 7. Leaves and Hodge side (s04, s06)

| Note | Macro | Meaning | LCY/RIG | LAP/LP2 |
|---|---|---|---|---|
| lam_j (bold) | `\lf_j(\psi)` | 1/(E - j psi), j >= 1; lam_0 = 1/E appears only in the lam_0 identity. In z = -psi, lam_j = 1/(E + jz) | lam_j | lam^E_j |
| L_m (script) | `\Lsp_m`, `\Lsp_m^{+}` | span(lam_1, ..., lam_m); constants + L_m | Lam_m | Lam_m |
| q_o | `\qo` | E - E^2 lam_1 = -E psi/(E - psi) = q_t at t = 0 | q_0 | q_0 |
| T_dil | `\Tdil(\psi)` | -psi^2/(E - psi) = q_o + psi | T | T^E |
| F^pot_h | `\Fpot_h` | genus-h potential of the Hodge CohFT Hg_w on the domain Dom | Fpot_h | Fpot_h |
| d^s (frak), Dom | `\dbase^s`, `\Dom` | base-point derivatives; q_o + t R[[psi]] | dd^s, Dom | 𝔡^s, 𝔇 |
| Xi_h(U) | `\Xiz_h(\msU)` | d_{lam_{j_1}} ... d_{lam_{j_r}} F^pot_h(q_o), scale E, class Hg_w | Xi_h(R) | Xi^E_h |
| Xi^tau_h | `\Xit_h(\msU)` | the scale-1 integral of [LP2] (quoted only) | | Xi_h |
| M, U (sans) | `\msM`, `\msU` | multiset of integers >= 2; of positive integers | R; U | R; U |
| dg, J | `\dg(\msU)`, `J_{\msU}` | Sum (j - 1); Sum j | D(R), J_R | D(U), J_U |
| Z[h, n] | `\ZZ[h,n]` | Xi_h(1^n); Z[h, 0] = F^pot_h(q_o) | Z[h, n] | Z[h,n] |
| M_irr(h) | `\Mirr(h)` | Xi_h(M) = 0 for every non-empty M with parts >= 2 and dg(M) <= 2h - 3 | M_irr(h) | Theorem L |
| tau_w, sigma_w | `\tauw`, `\st` | -epsilon_1/E; -epsilon_1 epsilon_2/E^2 (only when quoting [LP2]) | tau_w, sigma | tau, sigma |
| (H0), (S), (C), (L0) | text | dilaton homogeneity; string identity; vanishing on 1^p; the lam_0 identity | (H0), (S_P2), (C_P2) | LP2 s.7.1 |
| upsilon | `\upsilon` | the increment q - q_o in t R[[psi]] of a point q of Dom, in the definition of F^pot_h on Dom (def:set-fpot, s02) and in the Euler identity and the vertex-factor expansion of s06. Renamed on 2026-11-03: the sources write delta, which in the note is the algebraic function delta of def:set-alg | delta | delta |

## 8. Rigidity, three-torus, flop (s07)

| Note | Macro | Meaning | RIG |
|---|---|---|---|
| N (numerator) | `N` or `\mathsf{N}_{\sG}` | G_g = N/((1 - Q)^n Disc^m) in K_CY[Q]; if `N` clashes locally with the critical polynomial, use the subscripted form | N |
| ev | `\ev` | substitution Q -> Q_1 on K[[alpha_0]][[Q]] | ev |
| val | `\val` | alpha_0-adic valuation | val |
| N_d, L_d | `N_d`, `L_d` | genus-1 integrals of Lemma 5.7 | N_d, L_d |
| V_+-(h), Pi_+- | `V_\pm(h)`, `\Pi_\pm` | degree-1 localisation vertices; products of the five weights at p_+- | V_p(h), Pi_p |
| x (sans), c (sans) | `\xc`, `\cc` | alpha_1/alpha_2; 1 + x = E/alpha_2 on T_3 | x, c |
| kappa | `\kap` | epsilon_1 epsilon_2 alpha_2^4/alpha_1 (rescaling of the pairing on T_3) | kappa |
| Z', X', wt'_+- | `\Zfl`, `\Xfl`, `\wts'_\pm` | the flop and its tangent weights (alpha_1 - alpha_2, alpha_2, A, -epsilon_1, -epsilon_2), (alpha_2 - alpha_1, alpha_1, B, -epsilon_1, -epsilon_2) | Z', wt'_+- |
| s (local) | `\sQ` (= `\mathsf{s}`, main.tex) | 1/Q near Q = infinity (s02 lem:set-limits, s07, appendix; **not** the Euler variable, which does not occur there) | s |

## 9. Conventions where the note's reading differs from BS, or where care is needed (neutral wording)

State each once, where it first matters.

- d >= 0: F_g includes d = 0. BS's series has d > 0. Stated in s02 (definition of F_g); the difference is the constant F_g(0), which the rigidity statement subtracts.
- The conjecture of BS concerns beta != 0 and the sum over genera; the note proves the statement genus by genus for every d >= 1 and says in s07 exactly which group actions are covered (T_CY and everything acting through it), without listing what BS verified by computer except as a one-line remark.
- Genus 1 for d >= 2 uses a result of BS (the formula on the subtori alpha_i = epsilon_j). Say so plainly in s07; do not call it a gap.
- The localisation definition of the invariants (BS) is the one in force: the T_CY-fixed loci are the alpha_0-fixed loci of the proper moduli space of maps to P^1.
- Frobenius data at generic weights: the note derives them from the genus-0 potential (classical part by localisation plus the Aspinwall-Morrison term); BS's appendix is cited for comparison only, by content, with a \draftnote.
- Sign of z in the R-matrix and in the legs: as in the companion paper (legs carry R^{-1}(psi), rows v_i = e_i^* o R^{-1}(-z)); no comparison with other conventions beyond one sentence.

## 10. BS references: phrase to use, label for \draftnote

Line numbers never appear. Use the phrase, and `\draftnote{BS label: ...}` where the authors must check the pointer.

| Content | Phrase | BS label |
|---|---|---|
| definition of the equivariant invariants by localisation on the proper fixed locus; the fixed-degree and fixed-genus series | the definition of the equivariant Gromov-Witten invariants by localisation; the generating series | defn: refined GW invariant, defn: refined GW series |
| contravariance under group homomorphisms | the contravariance of the invariants | lem: GW contravariance |
| polynomiality epsilon_1 epsilon_2 GW in R_{G,2g} for proper moduli | the polynomiality of the invariants for a rigid class | lem: shape of refGW series |
| rigid class, equirigid target | rigid; equirigid | defn: rigid |
| the Rigidity Conjecture | the rigidity conjecture of [BS] | conj: rigidity |
| the unrefined case epsilon_+ = 0 of the conjecture | the case in which the action on X is also Calabi-Yau | prop: rigidity unref limit |
| closed formula on the subtori alpha_i = epsilon_j | the formula of [BS] for the resolved conifold on the subtori on which alpha_i = epsilon_j | prop:ResConifold, eq:resconGWref |
| the polynomials beta_n and their first values | the polynomials beta_n of [BS] | eq:betapol |
| the conifold asymptotics in the direct-integration axioms | the conifold asymptotics of [BS] | def:dirint, eq:bkgcf |
| the tangent weights at the fixed points, the Frobenius data of the conifold | the data of [BS, Appendix A] | App. A (no label) |
