# DIGEST_D2, part 5: s.3 Proofs, part B (formal stationary phase); s.4 Conventions

Continuation of `DIGEST_D2_4.md`. For each lemma of `DIGEST_D2_3.md` s.2B: the proof as numbered steps, where each step is written
out, and the recommended proof of record. Unless marked, every step is written out in full in the cited section.

## 3B. Proofs of the formal lemmas

**SP-DEF, Lemma 1.1** (FSL s.1.2). (1) Finiteness per weight bounds the weights below; products have finitely many contributing pairs per
monomial. (2) `d/dt_i`, multiplication by `t_i` act termwise. (3) `P^m` has weight `>= m/2`, so `exp(P)` has finitely many terms per weight.
`E_zC` preserves weight with finitely many inputs per output power (FSL s.1.3). RL s.1: PROVED.

**SP-U** (CP s.1.2 (F1); FSL s.1.5). Expand the definition: a term with `V` vertices of degrees `k_i` and an amplitude monomial of degree `j`
from `A_m` has `Sum k_i + j = 2(n - m + V)`, so `k_i <= 2n + 2`, `j <= 2(n - m)`; `[Phi_3^V]_d = 0` for `d < 3V`. (U1), (U2) are immediate.

**SP-H** (RL s.3.2). (1) `Hess(mu F) = mu H`, covariance `C/mu`, cubic part `mu F_3`. (2) For a monomial,
`E_(zC/mu)[t^alpha z^j] = mu^(-|alpha|/2) <t^alpha>_C z^(j + |alpha|/2) = <t^alpha>_C mu^j z'^(j + |alpha|/2) = E_(z'C)[t^alpha (mu z')^j]`. (3) Linearity and
weight preservation extend it to `R_n(A)`. (4) Take `P = e^(mu F_3/z) a(t; z) = e^(F_3/z') a(t; mu z')`. CP's two-line proof (CP s.1.2 (F3)):
the integrand depends on `(Phi, z)` through `Phi/z'`; `eta' = mu^(1/2) eta`.

**SP-IBP** (CP2 s.2.1; FSL2 s.3.1). (1) `T(eta) := exp(Phi_3(z^(1/2) eta)/z) h(p + z^(1/2) eta, z)`. (2) `(d_i Phi)(p + xi) = (H xi)_i + d_i Phi_3(xi)` and
`d/d eta_i = z^(1/2) d/d xi_i` give `exp(Phi_3/z) [(d_i Phi) h/z + d_i h] = z^(-1/2) [(H eta)_i T + dT/d eta_i]`. (3) Gaussian integration by parts
`<(H eta)_i P + d_i P>_C = 0` from the Wick recursion and `HC = -1`. (4) Apply to each coefficient of `z^(r/2)`. FSL's version: `g_i := e^(F_3/z) h_i`,
`d_i` a continuous derivation of `R_n(A)`, then (G3) (`E_zC[t_k g] = z Sum_j C_kj E_zC[d_j g]` by (W)); convergent sums by continuity.

**SP-PG.** *CP2 s.2.2 (purely algebraic; recommended):* (1) replace `delta` by `s delta`; `Lam_s(P)` := left side, `rho_s(P)` := right side;
(2) `d/ds Lam_s(P) = Lam_s(QP)`, `Q := (1/2) eta^T delta eta`; (3) `dC/ds = C delta C`, the double Wick step
`<QP> = (1/2) tr(delta C) <P> + (1/2) Sum (C delta C)_ij <d_i d_j P>`, `d/ds log d(s) = (1/2) tr(C(s) delta)` give `d/ds rho_s(P) = rho_s(QP)`;
(4) equal `s`-Taylor coefficients `<Q^k P>_C` by induction; `s = 1`; `I`-adic continuity. *FSL s.2 (G5):* both sides are polynomials over Q
in the entries of `C_0` coefficientwise in `k`; on real positive-definite `C_0` with `||k|| ||C_0|| < 1` they are Taylor coefficients of an
exact Gaussian integral (complete the square, dominated convergence), the only real integral of FSL s.1-s.8; Zariski density. (G5')
by rescaling and substitution `k = K_1`. Both PROVED (RF s.2; RL s.2 "the transfer is legitimate").

**SP-PROD** (CP2 s.2.3; FSL (G4)). Block-diagonal covariance; Gaussian moments of products of polynomials in disjoint blocks factorise.

**SP-IFT** (FSL2 s.3.2). `Psi(c) := c - Hess F(0)^-1 grad F(c)` maps `I^n` to `I^n`; `c - c' in I^r` implies `Psi(c) - Psi(c') in I^(r+1)`; iterates
converge, fixed point unique; `Hess F(c) = Hess F(0) mod I` is invertible.

**SP-T** (FSL2 s.3.2). (1) Admissibility of `b_p` in cases (a), (b). (2) `B(s) := b_(sp)`: each coefficient is a polynomial in `s` (Lemma U).
(3) `dB/ds = Sum_i p_i (d_i B + z^-1 (d_i F)(c + t) B)` (re-derived from the definition by RL s.1). (4) FL1 with `h_i = p_i B(s)` (coefficientwise
in `q` in case (a)): `d/ds SP^norm_c[F; B(s)] = 0`; a polynomial with zero derivative is constant.

**SP-T'** (FSL2 s.3.2). (1) `c` by the IFT over `A_0[[q]]`; `p = c - c_0`. (2) Lemma T(a) with amplitude `e^(delta/z) a`. (3) Expand both phases at
their critical points; the cubic parts `F_(0,3)` cancel, leaving `e^((F(c) - F_0(c_0))/z) e^(t^T (H_c - H_0) t/(2z)) e^(F_(c,3)/z)`. (4) (G5') with
`K_1 = H_c - H_0` (entries in `qA`). CP's (L1c)(iii) (CP2 s.2.6): iterate (iii').

**SP-PAR** (FSL2 s.4; CP2 s.2.6). (i) Lemma U, holomorphy of Taylor coefficients (Cauchy integrals) and of `C(lam)`; CP: `-H^-1 = -adj H/det H`
gives the rational form. (ii) (U1) for the Taylor map and for evaluation. (iii) Corollary T'. (iv) (1) Corollary T' with `q = eps`;
(2) every term of `eps g/z` has weight `>= 1/2` (`g = O(t)`), so `S(eps; z)` lies in the weighted ring with `eps` of weight 1 (RL-C4);
(3) the substitution `eps -> z` relabels coefficients, finitely many per weight, commutes with `E_zC`, and maps `e^(eps g/z)` to `e^g`.

**SP-DQ** (FSL2 s.4.5). (1) `q = q_0 + eps`; by (iii) at `q_0`, `S(q_0 + eps)` is a fixed functional of an `eps`-dependent amplitude; (2) `d/deps`
acts on the amplitude; (3) apply (iii) again. The iterated clause needs `D_q a = 0` (RL-C2). CP2 s.2.6 (iii'): translate to a fixed critical
point ((L1a), Jacobian 1); constant term, infinitesimal (PG) for the Hessian, cubic terms and amplitude differentiated; the `p'(t)`-terms
form `SP_(p_0)[Sum_i p'_i d_i(e^(Psi_t/z) B_t)] = 0` by (S1).

**SP-COV.** *FSL2 s.5:* (1) linear maps by (G2) (`C' = J^-1 C J^-T`); (2) `phi = J psi`, `D psi(0) = I`; (3) homotopy `psi_s = y + s chi`,
`M_s = I + s D chi`, `m_s := (b o psi_s) det M_s`, `Phi(s) := SP^norm_0[G o psi_s; m_s] = E_zC[e^(P_s/z) m_s]` with ONE fixed `C`, `P_s :=
G o psi_s - G(0) - y^T H y/2` (RL-C3); every coefficient is a polynomial in `s`; (4) with `W := M_s^-1 chi`: (a) `d_s(G o psi_s) = Sum_k d_k(G o psi_s) W_k`;
(b) `d_s m_s = Sum_k d_k(m_s W_k)` (Jacobi's formula and the Piola identity); (5) `Phi'(s)` is the SP of an exact form: 0 by FL1. Piola:
the `(n-1)`-forms `d psi_1 ^ .. ^d psi_i ^ .. ^ d psi_n` (factor `i` omitted) are closed; their coefficients are the minors. *CP2 s.2.4:* translations,
linear maps, `chi = M o chit`, the homotopy `chi_t(x) = x + t (chit(x) - x)` (Gaussian fixed, `d_t Phi_t = O(|x|^3)`), the vector field
`Y_t = (chit - id) o chi_t^-1`, Jacobi's formula and naturality of the Lie derivative turn `I'(t)` into the SP of an exact form; (S1).
Both PROVED (RF s.2; RL s.1, s.2). Either may be the proof of record.

**SP-CFU** (FSL2 s.5.4). L1a over `A = Q[[u]]` with `phi(t) := Y(1 + t; u) - Y_c(u)`; `phi'(0) = Y'(1; u) in 1 + u Q[[u]]` is a unit; the hatted
Hessian at `Y_c(u)` is a unit (`f''(1) = 2V(0) + 3 = 5` at `u = 0`); the critical values agree.

**SP-FUB.** *CP2 s.2.5:* (i) IFT; `v_c' = -H_vv^-1 H_vw`; `phi'' = S`; block determinant. (ii) (L1c)(i). (iii) (1) (L1a) with
`chi(v', w) = (v_c(w) + v', w)` (Jacobian 1): assume `v_c = 0`, so `H = diag(H_vv, S)`; (2) `Phi = phi(w) + (1/2) v^T H_vv(w) v + Phi_(3,v)`,
`delta(w) := H_vv(w) - H_vv(0)`; (3) the total expectation is iterated (block-diagonal covariance); (4) inner step: (PG) in
`R = C[eta_w][[z^(1/2)]]`, `I = (z^(1/2))`, gives `beta(w)` and the inner covariance `-H_vv(w)^-1`; what is left inside is the inner series
`B(w, z)` at `w = z^(1/2) eta_w`; (5) the outer expectation is `N[phi, beta B; 0]`; finitely many terms per order. *FSL3 s.6:* straighten
by `Theta(x', y) = (xi(y) + x', y)`; then (G4) and (G5') with `K_1 = K(y) - H_xx`; "Rings": the inner functional is identified with the
inner SP coefficientwise (RL s.2: `R_(n+m)(A)` is not contained in `R_n(A[[y]])`, e.g. `Sum_k x^(2k) y^(2k) z^(-k)`; RL-C5).
*Remark (FSL3 s.6):* no uniformity in the outer variables is needed; the inner data are power series in `y - y_c`.

**SP-L2.** *CP2 s.2.7 (degree count; recommended, needs no change of variables, RL s.2):* fix `w`; `s = sigma - sigma_c`, `xi = x_0 - X_0`,
`g_k := [d^k G/dx_0^k](X_0)/k!`, `fh_k := [d^k (sigma_c G + f)/dx_0^k](X_0)/k!`; `Phi - Phi(p) = s Sum_(k>=1) g_k xi^k + Sum_(k>=2) fh_k xi^k`;
`H = [[0, g_1], [g_1, 2 fh_2]]`, `det = -g_1^2`, `C = g_1^-2 [[2 fh_2, -g_1], [-g_1, 0]]`: `C_(xi xi) = 0`. Every monomial `s^a xi^b` of
`Phi_3^V A_m` has `a <= b/2`; a non-zero pairing needs `b <= a` (each `xi` pairs with an `s`); so `a = b = 0`: `N_m = A_m(X_0)`. *FSL3 s.7:*
change of variables `(s, g)` with `G(tau + T(g), w) = g` (L1a over `A_w`), `C_gg = 0`, every perfect matching of `g`'s contains `C_gg`;
item (4) by L1b, `(rho^2)^(-1/2) = rho^-1` in `1 + (w - w_c) A_w`.

**SP-FE** (FSL4 s.9.2). (1) `D = d/dY` raises the order in `1/Y`; power series in `D` act; `e^D f(Y) = f(Y + 1)`. (2) `D^2` is injective on
`Y^-1 Q[[1/Y]]` and `rho in Y^-2 Q[[1/Y]]`. (3) `D^2 mu = Sum_k B_2k Y^(-2k-1) = T - 1/Y - 1/(2Y^2)`, `T := (D/(e^D - 1))(1/Y)`. (4)
`(e^D - 1) D^2 mu = -1/(2 Y^2 (Y+1)^2)`. (5) `D^2[(Y + 1/2) log(1 + 1/Y)] = 1/(2 Y^2 (Y+1)^2)`. (6) `D^2 rho = 0`, so `rho = 0`. Bookkeeping for
`St_(m+1) = (1 + mw) St_m`: substitute `Y = 1/w + m` (FSL4 s.9.1).

**SP-L3.** *Proof A (FSL4 s.9.3; purely formal; proof of record per SYN5b K2):* (H) homogeneity `[z^k] SP^norm_a[F; (x/a)^e] = c_(k,e) a^-k`
(power count `r - |alpha| + |alpha|/2 = -k`), so `N_m = psi_m(z/a)`; A1 FL1 with `h = x^m` gives (v) and `N_(m+1) = (1 + mw) N_m`; A2 L1c(iv) with
`g = log(x/a)` and amplitude `(x/a)^-1`: `psi_0(w) = (1 + w)^(1/w + 1/2) e^-1 psi_0(w/(1 + w))`; A3 `Lam(Y) := log psi_0(1/Y)` satisfies
`Lam(Y) - Lam(Y + 1) = (Y + 1/2) log(1 + 1/Y) - 1`, the operator is injective on `Y^-1 Q[[1/Y]]`, `mu` solves the same equation (SP-FE), so
`psi_0 = St_0`; A4 all `m` (A1, SP-FE; `1 + mw` a unit); (ii), (iii) bookkeeping. RL s.1-s.2: complete; RL2 s.7.2 re-did A2 by a direct
Wick computation. *Proof B (classical; remark):* CP2 s.2.8 steps (1)-(5) (reduction to `a > 0` by (L1a) and (F3); Euler's integral;
Stirling DLMF 5.11.1; Laplace with the Morse substitution and Watson's lemma DLMF 2.3.8; uniqueness of asymptotic expansions) or FSL4
s.9.4 (DLMF (2.3.13)-(2.3.16), 5.11.1, s.2.1(iii)). RF-C6: the Pochhammer part is formal by (S1) with `h = x^M`; analysis needed only for
`M = 0`. *Conflict K2/O7:* CP2 s.2.8 "Not done: a purely formal proof" and CP4 s.5 item 4 are superseded by Proof A (SYN5b K2).

**SP-5118** (LF s.4.2 = LA s.2.3, same argument). (1) Both sides have coefficients in `Q[theta]`. (2) At `theta = 0` both are
`Sum_k b_k w^(1-2k)` (`B_(n+1)(0) = B_(n+1)`, `B_odd = 0` beyond 1). (3) Left side `Lhs(theta)`: with `log(w + theta) := log w + log(1 + theta/w)`,
`Lhs(theta + 1) - Lhs(theta) = log(1 + theta/w)` by SP-FE at `Y = w + theta` (continuous substitution `1/Y = w^-1 (1 + theta/w)^-1`). (4) Right
side: `B_k(theta + 1) - B_k(theta) = k theta^(k-1)` (DLMF 24.4.1; from the generating function 24.2.3, the difference being `t e^(theta t)`, LA s.2.3)
gives `Sum_n (-1)^(n+1) theta^n w^-n/n = log(1 + theta/w)`. (5) The difference `D_n(theta)` is a polynomial, periodic under `theta -> theta + 1`,
zero at 0: zero on `Z_(>=0)`, hence identically. RA's route (uniqueness of asymptotic expansions, DLMF 5.11.1 at a shifted argument)
is also valid (LF s.4.2); it is not formal.

**SP-L3C** (LA s.2.3). (1) By (F1), `[z^n] N[log t - t, t^(M-1); 1]` is a polynomial in `binom(M-1, j)`, `j <= 2n`, hence in `M`; so is
`[z^n] exp(Sum beta_n(M) z^n)`. (2) For `M = m in Z` the left side is `St_m` (SP-L3 (i)) and `log St_m(w) = Sum_n beta_n(m) w^n` is SP-5118 at
`w -> 1/w`, `theta = m`. (3) Polynomials agreeing on Z agree. (4) General `a`: CP2 s.2.8 step (1) uses no integrality. The RS critic
re-proved it by hand and VERIFIED it with a symbolic-`M` Wick sum (RS A3).

**SP-REFL** (RL s.3.3). SP-FE at `Y = -y`: `Gamma_St(1 - y) = (-y) Gamma_St(-y)`; sum the logarithms; `mu(-y) = -mu(y)`. The branch part
(`log y - log(-y) = i pi eps_j`, `e = 0`) is not a formal lemma (RL s.3.3; checked in RF s.1, RA s.1.1).

**SP-RES** (FSL4 s.10). Lemma 10.1 (closure of `L` under products: `ord_v >= min(eps', eps'')|j| - M' - M''`; `e^(-v f_0) in L`); FL2 (the
`X^-1` coefficient of a derivative is `0 . c_0`); FL3 (telescoping `Om(X + eps)(1 + eps') = Om(X) + d/dX [Sum_(j>=1) Om^((j-1)) eps^j/j!]`; the
Taylor shift is a continuous ring map, so the pull-back equals the naive one).

**Gaps in part B.** None in the lemmas the chain uses; the corrections RL-C1..C8 are statement precision only (RL2 s.5). Outside the
lemmas: the relabelling of FSL4 s.10.4 (SKETCH there; written in PB s.2.7, refereed once by PS; other digest). Not covered and not needed:
gamma's (R3) (FSL5 s.12).

## 4. Conventions, sign dictionaries, normalisations

* **q versus -q** (CP s.1.7; CP4 (B1); RA s.1.3; LA s.2.6). Every form has the exponentiated critical equation `-q (5 sigma)^5 =
  prod_j (zeta^j - sigma)`, i.e. `u sigma^5 = 1`; with `+q` the critical points would satisfy `sigma^5 (1 + 5^5 q) = 1`. `lambda = log(-q)` is forced
  by criticality. `c(z; L)` does not see `lambda` or the sign of `q`; they fix only the dictionary `L <-> q`. LMI's chart parameter
  `Q = e^(-lambda/alpha0)`: `Q = -q` at `s* = 1` (LF s.4.1), `Q = (-q)^(-1/5)` at `s* = 0` (LF s.4.2); a change of the sign of `q` multiplies `Q`
  by `e^(-i pi/alpha0)` (LA s.2.6). GEN's general summand: `Q = (-1)^(N + l_conc) z^(Sum c_f) q` (GEN s.2 Lemma 2).
* **Colours.** Colour `i` = critical point `sigma = zeta^i L` = fixed point `p_i`, `lambda_i = zeta^i` (CP s.0.2; LIT s.1.1). LP's convention:
  weights `lambda_0..lambda_4` on `H^0(P^4, O(1))`; tangent weights at `p_j` are `lambda_j - lambda_i`, `i != j` (LIT s.1.1).
* **Measures** (CP s.1.7; CP2 s.2.8; FSL4 s.9.1 (ii)). `Gamma(1 + w) = Int t^w e^-t dt` (measure `dy`, `M = 1`, for the twist variable);
  `Gamma(w) = Int x^w e^-x dx/x` (`M = 0`). Insertion `x^m` with `dx/x` is `M = m`, with `dx` it is `M = m + 1`. Literature (LIT3 s.4 item 6):
  Givental math/0108100 s.10 uses `dx`, Iritani Lemma 4.20 `dw/w`, GJR s.7 `dx/x` with `x^(k/m)` inserted. (Gamma) is blind to the
  twist measure (GAM-14); controls of the measure are meaningful only at `q != 0` (RF s.1).
* **Normalisations** (RL-C7; FSL5 s.12). `N` (= FSL's `SP^norm`) is Gaussian-normalised, `N_0 = A_0(p)`; `Nhat = N/N_0` has `c_0 = 1`. Symbols
  `(2 pi z)^(n/2) det(-H)^(-1/2) e^(Phi(p)/z)`; branch rules (mutually consistent, RL s.2): L1a `det(-J^T H J)^(1/2) := det J det(-H)^(1/2)`;
  L1b `Delta^(1/2) := Delta_in(y_c)^(1/2) det(-S)^(1/2)`; L2 (4) `(-G_t^2)^(-1/2) := (-g_c^2)^(-1/2) g_c/G_t`; Lemma H `det(-mu H)^(1/2) := mu^(n/2) det(-H)^(1/2)`.
  Only normalised series are canonical; `z`-dependent symbols must be explicit.
* **Phase conventions** (FSL s.1.4): `e^(Phi/z)` (CP, GAMMA) versus `e^(-W/w)` (POLY): `F = Phi`, resp. `F = -W`, `z = w`.
* **Stirling** (CP s.1.6; CP2 s.2.8): `log(a/z) := log a - log z`; `Gamma_St(1 + w) := w Gamma_St(w)`; `b_k = B_2k/(2k(2k-1))` (DLMF 5.11.1's
  coefficients); `St_0 = 1 + w/12 + w^2/288 - 139 w^3/51840 - 571 w^4/2488320 + 163879 w^5/209018880 + ..` (DLMF 5.11.4's `g_k`).
* **R-matrix and QRR** (LIT s.1.2-s.1.3; LIT3 s.4 item 6): Givental's `S ~ R e^(u/z) C`; arXiv:1712.02528 uses `R^-1` relative to
  Givental-Teleman (its footnote 11); math/0108100 Thm 9.1 (cotangent characters, sign +) agrees with math/0008067 Thm 1 (tangent characters,
  sign -); LP eq. (18) and ZZ eq. (5) are normalised by `Phi_s(0) = delta_(0,s)`, not by QRR; GJR's critical points `x_i = L - lambda_i`
  versus the chart's `x_j = zeta^j - L` ("up to `x -> -x`", LIT2 s.3.6). Twist rule (GEN s.3.1, s.4.1): an Euler-class twist contributes minus
  its fibre weight (convex `O(l)`: `-(l lambda_i + mu)`), an inverse Euler class plus its fibre weight (concave `O(-l)`: `mu - l lambda_i`).
* **Branches** (CP4 s.4 (B1)-(B10)): `lambda` forced; branch changes of `log y`, `log x_j` shift the phase by constants; `log x_0` is the
  equation `G = 0`; roots only in symbols; `e^(S_1)` normalised; `a^(M - 1/2)` only in symbols; `e = 0`; `L` independent (no fifth root);
  step (e) uses the power series `L(q)`. Conclusion: `c`, `Nhat_7`, `Nch` depend on `L` only.
* **Naming clashes to resolve in NOTATION.md:** the Wick functional `N` versus Conjecture N, N(g), the identity (N) and the numbers
  `N_(g,0)` (rename the functional); FSL's "(U)" (Lemma U) versus Theorem U; the DLMF shift `h` versus `h = h_(w_0)`; `Lam` (four meanings,
  D1) and `Lam1`; `Q` (LMI chart parameter) versus GEN's `Q` and the project's matrices `Q`, `Q^m`; `delta` (LMI shift, (PG) perturbation,
  T' perturbation); `K` (`= L - 1` in the chart) versus the slice `K`, `K_q`, a field `K`; `S` (LMI phase, Schur complement, fundamental
  solution); `C` (covariance) versus `C0`, `C_alpha`; `H` (Hessian, `H = d/dv`, hyperplane class); `c` (`c(z; L)`, Lemma N's `c(z)`, the
  source `c` of GAM-3); `g` (LMI phase, genus); `w` (hatted variable, weights `w_i`, Fubini's outer variables); `M0`, GEN's `M`, the CYC
  operator `M`; CP's (F1)-(F4), GEN's (L1)-(L4) and (C1)-(C4) versus correction tags; "(A2)" (GAM-7); PB's FL1-FL7 versus FSL's FL1-FL4.
