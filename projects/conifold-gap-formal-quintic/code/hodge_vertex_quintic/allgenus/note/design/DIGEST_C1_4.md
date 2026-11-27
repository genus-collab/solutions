# DIGEST_C1 (piece 4 of 6): s.3 Proofs

For each statement of pieces 2-3 that has a proof: the proof as numbered steps; after each step, where it is written out (file +
section). Marks: [W] written out in the source (possibly compressed: [W-c], expand in the note); [S] only sketched; [V] only
VERIFIED; [M] missing. Rule 2 requires every all-orders step to be written in sentences in the note; [W-c] steps need expansion,
[S]/[M] steps need a \draftnote unless another source supplies them.

## A. Tower and slice

**A1 Lemma 0 (frames).**
1. `<Y^C(-z), z^j phi_1> = z^j Y_1(-z) = -z^j + O(z^(j-5))` (from D6: `Y^C_1 = -1` plus powers `z^-n`, `n = 0 mod 5`, `n >= 5`). [W] RD1 s.1.3.
2. Existence for `z^k phi_a`, `a != 1`: `<Y^C(-z), z^k phi_a> = z^k Y_a(-z)` has no `z^k` term (`Y_a` has no `z^0` term for `a != 1`), so
   its non-negative powers are at most `z^(k-1)`; eliminate them from the top by subtracting multiples of `z^j phi_1`, `j < k`
   (step 1 makes each subtraction lower the top power). [W-c] RD1 s.1.3 ("eliminate non-negative powers from the top (`th_core.proj_T`)");
   the bound `j < k` is in the statement, the reason is implicit.
3. Uniqueness: a combination `Sum s_j z^j phi_1` in `K_0` has all `s_j = 0` (its top power survives by step 1). [W] RD1 s.1.3.
4. Direct sum `Hc_+ = K_0 (+) Q[[z]] phi_1` and "frame coordinates = flat components with `a != 1`": decompose each flat component with
   `a != 1` by step 2, keep the `phi_1` components. For a general element of `Hc_+` this is an infinite sum; at fixed `z^j` the
   corrections of `frame(k, a)`, `k > j`, carry `v^((4(k-j)+c_a+1)/5)` and converge u-adically (A3). [S] in RD1 s.1.3 (asserted with the
   statement); the convergence is in RD2 s.3.2 step (1) and TRF s.2. Writers: state the ring (homogeneous elements over `Q[v]`: finite
   sums; `Q[[u]]`: u-adic) and write this step.
5. Homogeneity of `f_j`: a weight count, `w(k,a) - w(j,1) = (4(k-j) + c_a + 1)/5 >= 1`. [W] TRF s.2 (Theorem R paragraph); S2A M5.

**A2 the slice.** 1. `Pi_i(-z phi_0) = 0` for `i >= 1`, `Pi_0(-z phi_0) = v`, `Pi_0(phi_1) = 1`, `Pi_i(phi_1) = 0` (`i >= 1`), from the
shape of `Y^C`. [W] RD1 s.1.3. 2. With A1: `{Pi_i = 0, i >= 1}` inside `-z phi_0 + Hc_+` is `-z phi_0 + (K_0 (+) Q phi_1)`, and `Pi_0` is
`v + s` on `-z phi_0 + n + s phi_1`. [W] RD1 s.1.3.

**A3 coefficient ring (M5).** 1. Exponent count of step A1.5. [W] S2A M5; TRF s.2. 2. Frame expansions converge u-adically and agree with
flat expansions inside a correlator (finitely many flat slots survive the dimension bound). [W-c] RD2 s.3.2 steps (1)-(2) parentheticals;
TRF s.2 ("the double sum converges u-adically and only finitely many flat slots are non-zero"). VERIFIED: `cr2_legdict` L4 (`k <= 12`).

**A4 model basis of K.** 1. `kappa(X, k) = [S(z) X z^k]_+` is the ancestor image of `tau_k(X)`, whose descendant times are `t^X_k`
only, so `Pi_i(kappa) = t^H_i = 0` for all `i`; hence `kappa in K_0`. [S] MS1 s.1 ("PROVED, elementary"; the descendant/ancestor
dictionary `[S^-1(z) q(z)]_+` is stated, not derived). 2. Log-freeness and non-negative powers: the columns `X != H` of `S` are log-free;
exponent `(c_a - c_(a+j) + 4j)/5 >= 1` for `j >= 1`. [W-c] MS1 s.1; MRF items 5, 8(c). 3. Triangularity `Pi_i(tau_k(H)) =
PiC(tau_(k-i)(H))`, diagonal `PiC(tau_0(H)) = 1`, gives `Hc_+ = K (+) span{tau_k(H)}`. [W-c] MS1 s.1. 4. Membership list from the support
rule `b + k = 1 mod 5` of `PiC` and the shift rule. [W-c] MS1 s.1. VERIFIED: MRF `rf_flat`.

**A5 consistency.** (a) `l_(-1) V^(k) = -V^(k+1)`: the `z^0` part of `q` drops out for `k >= 1` [W] RD1 s.2.1, TRF s.2; the string
counterexample [W] RD1 s.2.1, RD3 s.4.2. (b) cone and linearity [W] RD1 s.2.1. (c) `(2h-2) F_h + v F_h' = 0` for `F_h = c_h v^(2-2h)`; genus 1
`5/24 - 1/12` [W] RD1 s.2.1 (the values 5/24, -1/12 are computed: [V] MS1 s.2, TRM s.3). (d) [W] RD1 s.2.1, MS2 s.9(c). (e) [W] MS2 s.9(d).

**A6 isotropy.** 1. `Y^C(-z)^T eta^-1 Y^C(z)` is constant in `u` (flatness, `eta^-1` proportional to `eta`). 2. At `u = 0` it is
`e_1^T eta^-1 e_1 = 0`. [W] RD1 s.2.2; TRF s.2.

## B. Transport

**B1 Lemma E.**
1. Put `alpha(x) = eta f(-x)` (only negative powers of `x`); the contraction is `C(psi') = Res_x alpha(x)^T E_R(x, psi') dx`. [W] RD2 s.3.1.
2. Expand `1/(x + psi') = Sum_j (-psi')^j x^(-j-1)` in each of the two terms of the numerator; legitimate because their sum is divisible
   by `x + psi'`, which is (Sp). [W] RD2 s.3.1; TRF s.2 (says "legitimate in `Q((x^-1))[[psi1]]`"). Writers: fix the ring of the expansion
   and of the residue pairing (s.7 item O-8).
3. The `eta^-1` term: `Res alpha(x)^T x^(-j-1) = 0` (powers `<= -2`). [W] RD2 s.3.1; TRF s.2.
4. The other term gives `-Sum_j (-psi')^j [x^j](alpha(x)^T R^-1(x)) eta^-1 R^-1(psi')^T`, i.e. `C^T = -R^-1(psi') eta^-1 beta_(>=0)(-psi')`,
   `beta(x) = R^-1(x)^T alpha(x)`. [W] RD2 s.3.1.
5. Symplecticity: `R^-1(x)^T = eta R(-x) eta^-1`, so `beta(x) = eta (R f)(-x)` and `beta_(>=0)(-psi') = eta (R f)_(>=0)(psi')`; hence
   `C = -R(psi')^-1 (R f)_(>=0)(psi')`. [W] RD2 s.3.1; TRF s.2 ("signs, transposes and the substitution `x -> -psi1` are right").
6. "iff": `R(psi')^-1` is invertible. [W-c] implicit in RD2 s.3.1. Second slot by symmetry of `E_R`. [W] RD2 s.3.1.
7. Convergence: the infinite sums converge u-adically in the applications, `[z^-n] Y^C = O(v^(n - floor(n/5)))`. [W-c] RD2 s.3.1 remark.

**B2 Corollary E'.**
1. `E_P = E_R` for `P = P_0 R`: `P^-1(psi) eta^-1 P^-1(psi')^T = R^-1(psi) (P_0^T eta P_0)^-1 R^-1(psi')^T` and `P_0^T eta P_0 = eta` ((Sp)
   at `z^0`); `R` is symplectic for the same reason. [W] TRF s.2 ("because `P_0^T eta P_0 = eta`").
2. For `k >= 1`, `z^-k J^C in Hc_-` (`J^C` has a `z^0` term and no positive powers). [W] TRF s.2.
3. `(P z^-k J^C)_(>=0) = z^-k (P J^C)_(>=k) = 0` by (A); `(R f)_(>=0) = P_0^-1 (P f)_(>=0)` since `P_0` is constant in `z`. [W] RD2 s.3.1
   (first identity); the second is implicit [W-c].
4. Lemma E: the edges annihilate `Omega(z^-k J^C, .) = V^(k)` in each slot. [W] RD2 s.3.1.
5. Hence `E_P in K_0 (x) K_0` (completed): by A1 in each slot, a bivector annihilated by every `V^(k)` in a slot has no component along
   `Q[[z]] phi_1` there; write it in frames. [W-c] RD2 s.3.1 ("Hence") and s.3.2 step (1) ("write `E_P = Sum Et^(alpha beta) t_alpha (x)
   t_beta` with `t_alpha` the frames ... and `Et` the flat components with `a, b != 1`"). The duality step (the `V^(k)` restricted to
   `Q[[z]] phi_1` are triangular, so they detect the `phi_1` part) is implicit: write it.
Second proof: B5 (generating function). VERIFIED: s.5 EV-T1-EV-T4.

**B3 Lemma T.**
1. `Omega(f, P^-1 g) = Omega(P f, g)` by (Sp) (`Omega(P f, P g) = Omega(f, g)`); with `g = -z phi_0` this is (a). [W] RD2 s.3.1.
2. `f = z^-k J^C`: `V^(k)(q_T) = Omega(z^-k P J^C, -z phi_0) = -[z^(k-2)] Y'_0` (TRF form) `= (-1)^(k+1) [z^(k-2)] Y'_0(-z)` (RD2 form). [W]
   RD2 s.3.1; TRF s.2. Convergence of `Omega(f, q_T)` for `q_T in Hc_+` infinite in `z`: as B1.7 [W-c].
3. `k >= 3`: zero by (A); `k = 2`: `[z^0] Y'_0 = 0` because `Y'_0` has only powers `z^-n`, `n = 1 mod 5` (grading of `P` and D6); `k = 1`:
   `T_P`. [W] RD2 s.3.1; TRF s.2.
4. `q_T = -z P^-1(z) phi_0 = -(z/d) phi_0 + T_(>=2)`: in vector form `(P_0)^-1 phi_0 = phi_0/d` (`P_0` diagonal by the grading, `d = (P_0)_33 =
   (P_0^vec)_00`). [W-c] RD2 s.3.1 (the vector/lowered identification: TRF s.2, `P^vec = J P J`).
5. `qt = d q_T = -z phi_0 + d T_(>=2)` lies in `-z + K` (Lemma T(b) and linearity), `Pi_0(qt) = d T_P`. [W] RD2 s.3.1.
6. `s = O(u^2)`: at order `u^1 z^-1` in component 0 only `P_00(0, z) (-v'(0) e_0/z)` contributes, so `T_P = v'(0) (P_0)_00(0) u + O(u^2)`;
   `(P_0)_00 (P_0)_33 = 1` by (Sp) at `z^0`; so `d T_P = v'(0) u + O(u^2) = v + O(u^2)`. [W-c] RD2 s.3.1; TRF s.2 ("`[u^1 z^-1] Y'_0 =
   -v'(0) p_00(0)`"). The exclusion of the other `u^1` terms (e.g. `-u P'(0,z) e_1` has no `z^-1` term in component 0, by the grading) is
   implicit: write it.
7. Frame independence and the value `s = u^2/30 + (127/3150) u^3 + ..`: [V] TRF s.2-s.3 (`rf_flat`), RD2 s.3.3 (`gi_hat_table.log`).

**B4 Lemma F.**
1. Homogeneity at a vertex: `d^n Fb_h(q_T) = d^(2h-2+n) d^n Fb_h(qt)` (D4; for `h = 1` the `log` term is killed by `n >= 1` derivatives). [W]
   RD2 s.3.1 (the `h = 1` remark is implicit [W-c]).
2. Taylor expansion of `Fb_h` at the base point in the direction `d T_(>=2)`: each insertion has psi-degree `>= 2` and costs one dimension,
   so at most `3h - 3 + n` terms. [W] RD2 s.3.1; TRF s.2.
3. Only finitely many psi-degrees of the edge factors survive the dimension bound. [W-c] RD2 s.3.1.

**B5 model Lemma K.**
1. Contract `(w + z) V = eta^-1 - (R^-1(w) (x) R^-1(z)) eta^-1` in the first slot with `phi(w)`: `(w + z) C = eta^-1 phi(w) - R^-1(z) eta^-1
   ((R^-1(w))^T phi(w))`, whose right side has no positive power of `w`. [W] MS2 s.6.
2. Write `C = Sum_m C_m(z) w^m` (`C_m = Pi_m(V)` for `m >= 0`); the coefficient of `w^(m+1)`, `m >= 0`, gives `C_m + z C_(m+1) = 0`, so
   `C_m = (-z)^j C_(m+j)`; since `C_(m+j)` is a power series in `z`, the coefficient of `z^r` in `C_m` vanishes (take `j > r`). [W] MS2 s.6.
3. (b) is the definition of `Pi_m` applied to `-w R^-1(w) 1`. [W] MS2 s.6.
4. Application: `phi(w) <-> Y^C(-w)`; transpose form `(R^-1(w))^T = eta R(-w) eta^-1`; `Y^C_Q = -e_1 - t_c e_0/z + ..` gives
   `Pi_1(-z + T) = 0` and `Pi_0(-z + T) = t_c`. [W] MS2 s.6.
5. Convergence: in the flat frame `phi_j = O(u^(j - floor(j/5)))`; in the hatted frame the sums are finite at each `u`-order because hatted
   edges are polynomial in `(psi, psi')` ("Theorem D", not read). [W] MRF item 6 (correction of MS2's hypothesis).

**B6 Fock-space form.** 1. `Omega(f, R^-1 q) = Omega(R f, q)`. 2. By Lemma E, `exp(hbar E_R(d, d)/2)` commutes with multiplication by
`Omega(f, .)` (`[E(d,d)/2, L_f] = E(dL_f, d) = 0`). [W] RD2 s.3.1; TRF s.2. Quoted formula of arXiv:math/0108100 (no statement number;
conventions not re-checked).

**B7 the quintic case.** Inputs (G), Prop. S, Lemmas V, Q from GI s.2-s.4 (not read) and (U); then B1-B4. [W] S2A s.2 steps 1-4.

## C. Quintic Theorem R

**C1 Theorem R (a), (b).**
1. Edges: by B2 write `E_P = Sum Et^(alpha beta) t_alpha (x) t_beta` (frames `t_alpha`, `Et` = flat components with `a, b != 1`); the frame
   sums converge u-adically (A3). So every half-edge insertion is a frame. [W-c] RD2 s.3.2 step (1); TRF s.2.
2. Vertices of graphs with an edge: in a connected graph with an edge every vertex carries `n >= 1` half-edges, and every non-zero such
   graph has all vertex genera `<= g - 1` (a genus-g vertex forces a tree of genus-0 vertices, and a genus-0 vertex with `<= 2` half-edges
   vanishes by dimension because the leaves are `O(z^2)`). [W] RD2 s.3.2 step (2) (first sentence); TRF s.2 (the genus bound and the
   genus-0 vanishing); MRF item 8(a), (d).
3. Vertex factor (B3, B4): `d^(2h-2+n) Sum_m <t_(alpha_1) .. t_(alpha_n) (d T_(>=2))^m>_h / m!`, a finite sum. [W] RD2 s.3.2 step (2).
4. Decompose `d T_(>=2) = t' + s phi_1`, `t' in K_0` (frames with `k >= 2`, coefficients in `Q[[u]]`), `s in u^2 Q[[u]]` (A1, A2, B3). [W]
   RD2 s.3.2 step (2).
5. Multilinear expansion: every term is a model correlator with at least one `K_0` insertion and all others in `K_0 (+) Q phi_1`, times
   power series; by ML_K(h) it lies in `Q[v]` (the `F_h^(n) prod Pi_0` term vanishes because one `Pi_0` is zero). [W] RD2 s.3.2 step (2);
   MS2 s.7; MRF item 8(b).
6. Convergence of the frame expansion inside every correlator: frames with `k >= 2` have weight `>= 4/5`, so correlators with many of
   them carry high powers of `v`; the expansion agrees with the flat one (finitely many flat slots), and never splits a frame into its flat
   components, so the D-directions inside frames are never isolated. [W-c] RD2 s.3.2 step (2) parenthetical; MRF item 8(d) (in the flat
   frame each leaf of psi-degree `>= 2` raises `e` by `>= 3/5`; slot degrees at a vertex with `m` leaves `<= 3h - 3 + n - m`). MS2 s.10 calls
   these statements "argued, not formalised"; TRF s.2 checked them. Writers: write a lemma (u-adic convergence of frame expansions).
7. Edges are power series (flat regularity), the contraction is finite (B4): `Cont_Gamma in Q[[u]]`. This uses ML_K only for `h <= g - 1`.
   [W] RD2 s.3.2 step (2).
8. Vertex graph: `Cont_vertex = d^(2g-2) Sum_m <(t' + s phi_1)^m>_g / m!`; the terms containing `t'` lie in `Q[[u]]` by ML_K(g); the rest is
   `Sum_j s^j <phi_1^j>_g / j! = Sum_j s^j Phi_g^(j)(v)/j! + Q[[u]] = Phi_g(v + s) + Q[[u]]`, convergent since `s^j Phi_g^(j)(v) =
   O(u^(2-2g+j))`. [W] RD2 s.3.2 step (3).
9. Justification of the resummation: the sum over the number of leaves may be extended to all `m` (the added terms vanish by dimension)
   and re-ordered u-adically. [W-c] TRF s.2 (one sentence; not in RD2). Write it.
10. `v + s = d T_P` and `d^(2g-2) Phi_g(d T_P) = c_g T_P^(2-2g)`: (b). [W] RD2 s.3.2 step (3); TRF s.2.
Second proof (hatted, `P_can`): Givental's formula for `P^ A^LG`; Lemma K; ML in function form; `Pi_0` is linear and vanishes on `K_0`, so
`F_h(Pi_0(-z + T + x)) = F_h(t_c)` for `x in K_0` and every derivative in an edge direction only sees `Reg_h`; u-adic convergence; vertex
graph `F_g(t_c) + Reg_g(-z + T)`. [W-c] MS2 s.7; MRF item 8 (a)-(d). (c): one line [S] RD2 s.3.2 step (4).

**C2** is a reading of the proof C1 [W] RD2 s.3.2.
**C3 Cor. L.** 1. First-order Taylor term `c_g (2 - 2g) vhat^(1-2g) (t_c - vhat)`. [W] RD2 s.3.3; TRF s.2. 2. `t_c - vhat = u^2/30 + O(u^3)`
[V] only (`gi_hat_table.log`; TRF `rf_flat`; frame independence by `cr_conv` (4), SYN01 C6). No written derivation in my sources [M].
**C4 Prop. I.** Apply C1(c) to `P_(A') P_A^-1` (`P_A^-1 Y^C_A = Y^C`). [W-c] model-outwards; [S] otherwise. RD3 s.4.3.
**C5 deformation route.** (D0)-(D3) by differentiating `z P_s' = A_s P_s - P_s At`, Prop. S and Theorem G; `B_s` is the `epsilon`-derivative
of `E_R` at `R = exp(epsilon X_s)`. [W] RD3 s.4.1; TRF s.2. The Hessian argument [W] RD3 s.4.2.

## MOD, TR, UNI

**MOD-1** on K: the `n = 1` case of ML_K; a singular value with zero `Q[v]` part is the law [W] (follows from S2A s.1; MS1 s.3 "the one-point
law is exactly (ML) restricted to one point"). Off K: [V] exact (MS1 s.3; MRF item 2).
**MOD-2** [V] exact counterexamples (MS1 s.4; MRF item 3); the logic "on singular pairs the Hessian must be genus independent" [W] TRF s.4.
**MOD-3** 1. quote arXiv:2211.12259 Thm 2; 2. dictionary [W] MS2 s.8, with MRF item 5 (the constant `c = 0`); 3. Euler's formula count [W]
MS2 s.8; MRF item 5 (checked).
**MOD-4** [W] MS2 s.9(b) (formal; the Barnes difference equation is classical, not proved in the source).
**MOD-5** [W] MS2 s.9(c); MRF item 9.
**MOD-6** formula: first order of `exp(c z^m)` [W] MRF item 10 (re-derived); values [V] MS2 s.9, MRF item 10.
**TR-1** 1. hypotheses of arXiv:1509.06954 Thm 8.1 (critical points `p^5 = v/4`, `lambda'' = 20 p_c^2`, distinct critical values `5 p_c^4`)
[W] TRM s.1; TRRF claim 1. 2. TFT normalisation: `Delta_j^TR = Lambda c_j^2`, `e_j = sqrt(Lambda) et_j`, TR TFT `= Lambda^(g-1)` frozen;
`xi^j = xit^j/sqrt(Lambda)`, `xit^j = -1/(c_j (p - c_j))` [W] TRM s.2.2; TRRF s.3.1. 3. a CohFT rescaled by `Lambda^(g-1)` with metric
`/Lambda` is a CohFT (vertices `Lambda^(g_v-1)`, edges `Lambda`, total `Lambda^(g-1)`); the R-action commutes [W] TRRF s.3.1. 4. `R_TR = M_0(z/kappa)`
multiplies the degree-D part by `kappa^-D` [W] TRRF s.3.1. 5. units [W] TRRF s.3.1. 6. TR-2 [W-c]/[V]. 7. `n = 0`: TR-3 [W].
**TR-2** 1. (2-31) computed: [V] `z^0..z^8` (TRM `tr_run_rmat`), `z^10` (TRRF `rf_t5_rmat`). 2. Uniqueness of the homogeneous R-matrix:
[S] TRM s.2.3 (two sentences; "standard, not re-proved", TRM s.8); inputs to be named (TRRF correction 4). Writers: write the recursion
argument or cite it with a statement number; else \draftnote. Alternative direct route [S] TRRF correction 4.
**TR-3** integration by parts at each branch point; `P_(i,d) = Res_0 + Res_inf` of `xit^i_d y dx`, `y dx = kappa (p^4 - 1/p) dp`; at `p = 0`:
`-kappa xit^i_d(0)` (`1/x'` vanishes at 0 for `d >= 1`); at infinity `xit^i_d = O(p^(-1-4d))`. [W] TRRF s.3.3.
**TR-4** route 1: (i) quote (1.3), b-points, hypotheses [W] TRM s.4.1, TRRF s.3.2; (ii) extraction `I_g` [W] TRM s.4.1 (re-derived TRRF);
(iii) infinity: degree `<= 1 - (2g-2) r - 2g <= -5` [W]; (iv) the pole: `D^n(z) = (-1)^n n! z^(n+1) v^-n (1 + O(z^2))`,
`x^(2k) = (2k)! v z^-(2k+1) (1 + O(z^(2k+1)))`, lowest order exactly `z^-1`, `Res = v^(2-2g) A_g` [W]; (v) binomial resummation,
`1 + phi = artanh(s)/s`, `t = 4 s^2`, Lagrange inversion `[s^1] w^(1-2g) = (1-2g) [w^2g](w coth w) = (1-2g) 4^g B_2g/(2g)!`,
`A_g = -B_2g/(2g)` [W] TRM s.4.1; TRRF s.3.2. Route 2: TR-8. Independent short proof [W] XYRF s.2.
**TR-5** [W] TRM s.4.3 (with the "shifted family" input [S]: standard, not re-proved).
**TR-7** structure of (7.9): Laurent polynomials over powers of `4 z_i^5 - v` and `z_i - z_j`; no diagonal poles; regular at `z_i = 0` for
`v != 0` [W-c] XY s.2.1; XYRF s.3 ("`z_i`, `z_i - z_j` are monic").
**TR-8** [W-c] XY s.2.3 (monomial count, Borel sum in one paragraph); contour lemma [W] XYRF s.2 (two lines). The Borel step is a formal
identity of coefficient sequences justified by a substitution under an integral: write it as such, or present route 1 only.
**TR-9** 1. for legs in `L` the pairing is `(-1)^n Res_inf .. Res_inf prod Phi_j omega_(g,n)` (no poles at 0; functions of `x` drop out) [W-c]
XY s.3; 2. `omega_(g,n) = N/prod q(z_i)^K` with `N` polynomial in `z` and the parameters, and `1/q(z)` expands at infinity with polynomial
coefficients [W] XYRF s.3. 3. hence a polynomial in the parameters [W].
**UNI-1** [W] RD2 s.3.3, TRF s.2 (inputs (R1), Prop. S for `r != rt` from UR/GI [not read]); then C1 [W].
**UNI-2** (a)-(e) [W] BVRF s.3 (chain rule; integral substitution; pointwise construction; `ord_0` preserved); hatted count [W] BVRF s.3.
**UNI-4** [W] UN1 s.2.3; URF s.2.1. **UNI-5** [W] UN1 s.2.3; URF s.2.4. **UNI-6** [W] URF s.3.4. **UNI-9** [W-c] UN2 s.3.1, s.3.5.
**UNI-3, UNI-7, UNI-8**: statements about examples: [V] only (s.5).
