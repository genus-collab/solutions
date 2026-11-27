# DIGEST_C2, part 5: s.3 Proofs, continued (log legs, ML_K, S15, the chain, universality)

Continues `DIGEST_C2_4.md`; same markers [W], [W-], [Q], [V], [S], [M].

## `lem:mlk-fixedx` [all W; LL s.2.1; LR s.3 re-derived]
1. At fixed `x`: `d_v p|_x = -(d_v x)/(d_p x) = -(1/p)/((p^5 - v)/p^2) = -p/(p^5 - v)`.
2. For any family of forms, `d_v|_x = d_v|_p - Lie_xi` in each variable.
3. `d_v(y dx)|_x = (d_v p|_x) dx = -dp/p`; and `Int_{q in gamma} B(p, q) = dp [1/(p - q)]_{q=0}^{q=oo} = -dp/p`. So `d_v(y dx)|_x = Int_gamma B(., q)`.

## `lem:mlk-ibp` [all W; LL s.2.2; LR s.3]
1. (a) On the closed circle `Lie_xi(Phi w) = d(i_xi(Phi w))`, so `oint Phi Lie_xi w = -oint (xi . dPhi) w`.
2. Divide `p Phi'(p) = Q (p^5 - v) + r` (`deg r < 5`, monic divisor): `xi . dPhi = Q + r/(p^5 - v)`; `oint (r/(p^5 - v)) w = 0` since the
   integrand is `O(dp/p^3)` at oo; so `oint Phi Lie_xi w = -<(xi dPhi)_+, w>` (the constant term of `Q` pairs to zero).
3. (b) `Int_gamma Lie_xi w = [i_xi w]_0^oo = 0`: `xi(0) = 0`, `w` holomorphic at 0; at oo `xi = O(p^-4)`, `w = O(dp/p^2)`.

## `thm:mlk-var` ((Q1) as used)
1. EO Thm 5.1: (5-14) `D_Omega W^(g)_k(p_1..p_k)|_(x(p_i)) = Int_dOmega W^(g)_(k+1)(.., q) Lambda(q)` for `g + k > 1`; (5-15)
   `D_Omega F^(g) = -Int_dOmega W^(g)_1 Lambda` for `g >= 1`; temperatures, s.5.3: `Omega = -dS = Int_alpha^alpha' B`, `Lambda = 1`,
   (5-27)-(5-30) [Q; [DW] confirmed on the printed pp. 35-36].
2. Dictionary [W; LR s.2]: engine kernel = minus the printed (4-3) kernel, so `omega_{g,n} = (-1)^n W^(g)_n`, `F_g = -F^(g)_EO`; EO (5-1)
   `delta_Omega y|_x dx = -Omega` makes the engine hypothesis `d_v(y dx)|_x = +Int_gamma B` read `Omega = -Int_gamma B`, `Lambda_EO = -1`;
   then `d_v omega_n = (-1)^n (-1) (-1)^(n+1) Int omega_(n+1) = +Int omega_(n+1)` and `d_v F_g = +Int omega_(g,1)`.
3. Hypotheses (list in part 2) [W; LL s.2.3, LR s.2]; LR checked the Rauch case (5-4) independently. Slip in LR s.2: "`Omega = -dp/p`" in
   the hypothesis list, `Omega = +dp/p` in the Rauch check; by (5-1) `Omega = -d_v(y dx)|_x = +dp/p` [DW]; regularity of `Omega/dy` is
   unaffected.
4. Case `(0,1)`: the hypothesis itself; case `(0,2)`: `lem:mlk-genus0` [W].

## `lem:mlk-genus0` [W, terse in places; LL s.2.4; LR s.3 re-derived]
1. `d_v|_p B = 0`, so `d_v|_x B = -(Lie_xi1 + Lie_xi2) B = -d1 d2 [(xi(p1) - xi(p2))/(p1 - p2)]` (`xi` read as `p/(p^5 - v)`) [W].
2. Partial fractions `xi(p) = Sum_b r_b/(p - c_b)`, `r_b = 1/(5 c_b^3)`, no polynomial part; `(xi(p1) - xi(p2))/(p1 - p2) = -Sum_b r_b/((p1 - c_b)(p2 - c_b))`;
   hence `d_v|_x B = Sum_b r_b e_{b,1}(p1) e_{b,1}(p2)` [W].
3. `omega_{0,3} = -Sum_b e_{b,1}^3/x''(c_b)` from the leading order of the kernel (`K = e_{b,1}(p0)/(2 dx)(1 + O(t))`, `B(sigma q, .) = -B(q, .)(1 + O(t))`)
   [W-; one sentence in LL s.2.4].
4. `Int_gamma e_{b,1} = -1/c_b`, `x''(c_b) = 5 c_b^2`; so `Int_gamma omega_{0,3}(p1, p2, .) = d_v|_x B` [W].
5. The four (0,3) pairings by residues at the `c_b` and `Sum_b c_b^N = 5 v^(N/5) [5|N]`; for `Phi_i in p Q[v][p]` the powers of `c_b` are
   `>= -4`, only multiples of 5 survive, so every pairing with an L-leg is in `Q[v]`; `<(-log p)^3>_0 = (1/5) Sum c_b^-5 = 1/v` [W].

## `prop:mlk-V1` [all W; LL2 s.3 proof; LR s.3]
1. Fix `v_0 != 0`; for `v` near `v_0` choose `R > max_b |c_b(v)|` and a ray `gamma` avoiding all `c_b(v)`; then
   `C(v) = (2 pi i)^-a oint_{|p_i|=R} Int_{gamma^k} prod_i Phi_i(p_i) omega_{h,n}` with contours independent of `v`.
2. Differentiate under the integrals: holomorphic in `v` on the compact contours; on the rays bounded uniformly near `v_0` by
   `const |dq|/(1 + |q|)^2` (by (F1) and continuity of the coefficients). So `d_v C = Sum_i <.., d_v Phi_i, ..> + Int prod Phi_i d_v|_p omega_{h,n}`.
3. `d_v|_p omega_{h,n} = Int_gamma omega_{h,n+1}(., q) + Sum_l Lie_xi_l omega_{h,n}` (`lem:mlk-fixedx`, `thm:mlk-var`).
4. The first term is the pairing with one more log leg (Fubini; no diagonal poles by (F1)); a Lie term in an L-slot gives
   `-<.., (xi dPhi_i)_+, ..>` (`lem:mlk-ibp`(a)); in a log slot it gives 0 (`lem:mlk-ibp`(b)). Rearranging gives (V1).

## `thm:mlk-logleg` (Theorem V)
1. (i) Fix `h` and `a >= 1`; induction on `k`. Base `k = 0`: `thm:mlk-lreg` (for `2h - 2 + a > 0`) [W; LL2 s.3].
2. (i) Step `k -> k+1`: if `(h, a + k) = (0, 2)` the new pairing is a (0,3) pairing (`(a, k+1) = (2,1)` or `(1,2)`), explicit and in `Q[v]`
   (`lem:mlk-genus0`); otherwise every term on the right of (V1) is a pairing with `a >= 1` L-legs (`d_v Phi_i` and `(xi dPhi_i)_+` lie in
   `Q[v][p]` since `p^5 - v` is monic) and `k` log legs, or the `v`-derivative of one: in `Q[v]` by induction [W; LL2 s.3]. The bases for
   `(h,a) = (0,1), (0,2)` are the explicit (0,3) pairings [W-: LL2 phrases them as the step from the unstable level].
3. (ii) For `a = 0`, (V1) reads `<(-log p)^(k+1)>_h = d_v <(-log p)^k>_h` [W].
4. (ii) Base `h >= 2`: `<-log p>_h = d_v F_h^TR` -- by (Q1) at `n = 0` (LL2 s.3), or without (Q1) (LR F2): the definition gives
   `(2-2h) F_h = (1/5)<p^5>_h + v <-log p>_h` (`Phi = p^5/5 - v log p`), and homogeneity gives `<-log p>_h - d_v F_h = -(1/(5v)) <p^5>_h`, which is 0
   because `<p^5>_h in Q[v]` (L-regularity) has exponent `2 - 2h < 0`; then `F_h^TR` by `thm:mlk-gap` [W].
5. (ii) Base `h = 1`: homogeneity gives `<-log p>_1 = c/v` for one rational `c`; `c = -1/12` computed exactly over `Q(zeta5)` [V; LL2 s.3,
   `ll_falsifier_exact.log`]. Written alternative in the sources: `lem:mlk-vanper` at `g = 1` (`-(B_2/2) v^-1 = -1/(12 v)`), PROVED mod
   arXiv:2212.00320 Thm 7.1, applied to the form D1 [Q + W; XY s.2.3, XYR s.2; the transfer between forms: s.7 item W5].
6. (ii) Base `h = 0`: `<(-log p)^3>_0 = 1/v = F_0'''` for `F_0 = (1/2) v^2 log v` (`lem:mlk-genus0`) [W].
7. (iii) is `prop:mlk-V1` [W].

## `thm:mlk-main` (ML_K) [W; SYN2 s.2 step 7.1-7.6; LL2 s.3 Corollary with LR F3]
1. (M1) `<q_1 .. q_n>_h = (-5)^(1-h) <Phi_1 .. Phi_n>^TR_h` whenever `q_i = q_{Phi_i}`.
2. (LD) The TR leg of `s in S_w` is `-Lambda(s)`.
3. (D(ii)) The log coefficient of `-Lambda(s)` is `-Sum_i Pi_i(s)(-x)^i/i!`.
4. (M4, M5) For `X in K cap S_w`: `X = Y + Pi_0(X) tau_0(H)`, `Y in K_0 cap S_w` (Cor. 5); `leg(X) = -Lambda(Y) + Pi_0(X)(-log p)`; by D(iii)(b)
   and Lemma 1(a), `-Lambda(Y) = Phi + f(x)`, `Phi in Q[v][p]`, `f in Q[v][x]`.
5. (M6, M7, M9(i)) Multilinearity: `<X_1..X_n>_h = (-5)^(1-h) Sum_S prod_(i not in S) Pi_0(X_i) <(Phi_i)_(i in S), (-log p)^(n-|S|)>^TR_h`;
   every term with `S` non-empty lies in `Q[v]` (the factors `Pi_0(X_i) in Q v^(w_i/5)` have non-negative exponents).
6. (M8, M9(ii)) The term `S = {}` is `(-5)^(1-h) (F_h^TR)^(n) prod Pi_0(X_i) = F_h^(n) prod Pi_0(X_i)` (frozen `F_0 = -(5/2) v^2 log v`,
   `F_1 = -(1/12) log v` included); for `n = 0`, `F_h = c_h v^(2-2h)`.
* Note: "homogeneous with coefficients in `Q[v]`" means `X_i in K cap S_(w_i)`; insertions of quintic Theorem R are `Q[[u]]`-combinations of
  these (`prop:trm-coeff`), and the statement is applied termwise [W-; SYN2 M5].

## S15 lemmas and Theorem S15 [W; S15 s.2-s.4; FR s.1 line by line]
* `lem:flh-frame`: insert `A A^-1` on every half-edge and leaf slot (every term is a full contraction) [W].
* `lem:flh-dil`: `F^(kappa.Omega)_h(q) = kappa^(2h-2) F^Omega_h(kappa q)` termwise; so `d^n F^(kappa.Omega)_h(q) = kappa^(2h-2+n)(d^n F^Omega_h)(kappa q)
  = (d^n F^Omega_h)(q)` by the dilaton identity [W; FR adds the reading "q absolute"].
* `lem:flh-hat`: (i) `S = diag(s^(1-c_a))` has `s^(1-c_a) s^(1-c_b) = 1` for `a + b = 3 mod 5` (`c_a + c_b = 2`), so edges are unchanged; the unit
  row `a = 0` gets `s`, so every leaf gets `s`; (ii) with `z = lambda zhat`, `lambda = s^4`: vertex classes of degree `D` give `lambda^D`, edge
  coefficients `psi^m psi'^n` give `lambda^(m+n+1)`, leaf coefficients `z^k` give `lambda^(k-1)`; the degree at a vertex is `3h_v - 3 + n_v + j_v`
  (`j_v` leaves); total `lambda^(Sum_v (3h_v - 3 + n_v) + E) = lambda^(3g-3)`; (iii) `xhat^(2h-2+n+j)` gives `s^-(2h-2+n+j)` per vertex, the leaf
  factors cancel `s^-j`, leaving `s^-(2g-2)`; total `s^(12(g-1) - 2(g-1)) = u^(2g-2)` [W; applies to the hatted sum as implemented: s.7 item G4].
* `lem:flh-comm`: the inner Givental sum of `R.(kappa.omega)` is `Cont[E_R, R^-1(-z kappa^-1 1); kappa.omega]` with legs; each inner vertex gives
  `kappa^(2h_v - 2 + n_v)` (legs included), `Sum_v (2h_v - 2 + n_v) = 2g - 2 + n` [W; FR: the translation for `kappa.w` is `kappa^-1 T`, PPZ Def. 2.5].
* (I1)-(I3): from UR s.2, s.4 and GI Prop. S [not digested here; UR not read]. (I4): units [W; S15 s.3].
* `thm:flh-main`:
  1. By (F), `Q' = Q M_nat^-1 = P Q^m(v; z) M_nat^-1 = P Psi^m`, i.e. `R' = P o Psi^m` [W].
  2. `Psi^m` is an isometry: `E_{R'} = (Psi^m^-1 (x) Psi^m^-1) E_P`, `q_{R'} = Psi^m^-1 q_P` [W].
  3. `T_Gamma = Cont[E_R', q_R'; M_nat.omega^Q]` (1.3) `= Cont[E_R', q_R'; M_nat.(d.omega^m)]` (I4) `= Cont[E_R', q_R'; d.(M_nat.omega^m)]`
     (Lemma 4, `M_nat(0) = 1`) `= Cont[E_R', q_R'; M_nat.omega^m]` (Lemma 2) `= Cont[E_P, q_P; Psi^m_*(M_nat.omega^m)]` (Lemma 1)
     `= Cont[E_P, q_P; Omega^mod_{v(u)}] = F_Gamma` (I1) [W].
  4. `G_Gamma = u^(2g-2) T_Gamma` (Lemma 3) [W]. The unit hypothesis of D12 holds with `kappa = 1` for `R'` and `kappa = d` for `P` (`P_0`
     diagonal by the grading) [W].
  5. Window: under (F_g) instead of (F), step 1 holds modulo `z^(3g-2)`, which is all that a genus-g graph sees [W-; FR (I2); SYN2 s.2].
  6. u-adic convergence of the dilaton re-expansions: `d(0) = 1` [W-; FR s.1 (1.1 row), C5: "unstated"].

## The chain (steps 1-5: Stage-1 sources, not digested; steps 6-8 above and below)
1. Steps 1-3 [W in GI s.2-s.4]; step 4 [W in RD_2 s.3.1]; step 5 [W in RD_2 s.3.2 and MSG_2 s.7; "refereed twice"] (SYN2 s.2).
2. Step 8: combine step 5 (`F_Gamma in Q[[u]]` for edge graphs; vertex `c_g t_c^(2-2g) + Q[[u]]`, given ML_K(h) for `h <= g`), step 7 (ML_K for all
   `h`) and step 6 (`G_Gamma = u^(2g-2) F_Gamma` under (F_g)): every edge graph of the hatted outer sum is regular and the vertex graph is
   `c_g t_c^(2-2g)` + regular [W; SYN2 s.2 step 8].
3. The (F_g) input is discharged in the main theorem by the Stage 4-5 chain (MS5 s.2 items 1-6) [outside this digest].

## Universality
* `lem:uni-pole`: (i) the `z^0` equation of the gauge recursion with `P_(-1) = 0` and `P_0 = diag(q^0)` invertible [W; BR s.1(i)];
  (ii) `ord_0(v'^5/v) = 5(n - 1) - n = 4n - 5` [W; BV s.1, BR s.1(ii)]; `v'(0)^4 = -res_0 B` [W; BR s.1]. Converse paragraph [S; not needed].
* `lem:uni-coord` [W; BR s.3]: (a) `prod b~ = phi'^5 B_X(phi) = -phi'^5 v_X'(phi)^5/v_X(phi) = -(v_X o phi)'^5/(v_X o phi) = -v_q'^5/v_q = B_q`, and
  `b~_a = phi' b_a(phi)` keeps the shape; (b) chain rule `z d/dw P(phi) = phi'(N P - P Nt)(phi) = N~ P(phi) - P(phi) Nt~`, normalisations are
  conditions at the origin (`phi(0) = 0`); (c) `y_(n+1) = Int_0 N y_n` and `Int_0^(phi(w)) h(u) du = Int_0^w h(phi) phi' dw`; (d) `M = M_0(z v^(-4/5))`
  depends only on the model coordinate and `v_X~ = v_X o phi`; `R = Psi^-1 P Psi_mod R_mod(v)` pointwise; an amplitude is a polynomial in the
  coefficients of `R'`, `M`, the leaves and the TFT at the same base point; (e) `phi(0) = 0`, `phi'(0) != 0` preserve `ord_0`. Hatted encoding:
  `lam := (w/phi(w))^(1/5)`; rows of `X~` at `w` = `lam^(1-c_a+4n)` times those of `X` at `phi(w)` (`[zhat^n]`); `E_M` likewise
  (`tau_X~(w)/tau_X(phi(w)) = lam^4`); an edge `(m,n)` gives `lam^(4(m+n+1))`, a leaf `T(e)` gives `lam^(4e-3)`, a vertex `lam^-(2g_v-2+n_v)`;
  with `Sum psi-exponents = 3g - 3 - E + L` the total is `lam^(10(g-1))` [W, engine-level].
* `cor:uni-flat`: steps 1-5 and 7 for the system (Theorem U for the system) [W-; SYN2 s.2 Corollary, five lines]; Lemma C for B-perturbed
  systems [W]; the hatted-engine form via step 6 [S: "matching the univ engine's conventions to Lemma 3 has not been written out"].
