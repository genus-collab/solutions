# DIGEST_C2, part 4: s.3 Proofs, model side (TR identification, LD, D(ii), D(iii), L-regularity, Gap)

Continues `DIGEST_C2_3.md`. Each proof is a numbered list of steps; every step names where it is written out. Markers:
**[W]** written out in the source; **[W-]** written but terse (one clause; the writer must expand, transcribing only);
**[Q]** a quoted theorem; **[V]** only VERIFIED (finite computation: may not appear in a proof, binding rule 3);
**[S]** only sketched; **[M]** missing. Items marked [V], [S], [M] are collected in s.7.

## `thm:trm-ident` (M1)
1. DNOPS Thm 8.1 applies at every point `a = (0,0,0,0,v)`, `v != 0`, of the Hurwitz space `H_{0;1,4}` (hypotheses: `dig2` entry) [Q; TRM s.1].
2. (2-35): `omega_{g,k} = Sum Int alpha_{g,k}(e_i1..e_ik) prod psi^(d_j) d((d/dx)^(d_j) xi^(i_j))`, `e_i` normalised canonical [Q; TRM s.1].
3. Dictionary to the frozen theory: (2-30) `Delta_j^TR = Lambda c_j^2`, `Lambda = -5 kappa`; frozen `Deltahat_j^(1/2) = zeta^(4j) = c_j`, so
   `e_j = sqrt(Lambda) et_j` and the TR TFT on the `et_j` is `Lambda^(g-1)` times the frozen TFT; (2-34) `xi^j = xit^j/sqrt(Lambda)`,
   `xit^j(p) = -1/(c_j (p - c_j))`; flat basis `H^a = N_a Sum_i zeta^(i(a+1)) et_i`, `N = (1,-1,-1,1,1)`; `xif_a` = xi-function of the
   dual vector `H^(3-a)/5` [W; TRM s.2.2, "PROVED, elementary"].
4. R-matrix: (2-31) gives `X(z) = M_0(z/kappa)^-1 = M_0(-z/kappa)^T`, i.e. `R_TR(z) = M_0(z/kappa)` [V through `z^8`, TRM s.2.3]. All orders:
   both are homogeneous R-matrices of the same `(U, V)`; `[R_(m+1), U] = (m + V) R_m` determines `R` from `R_0 = 1`; `U` agrees (critical
   values); `V = [R_1, U]` agrees "because `R_1` agrees" [W for the uniqueness argument; the agreement of `R_1` is [V] (the `z^8` check);
   TRM s.2.3 "PROVED modulo that standard uniqueness"]. The sign of `kappa` (stationary phase of `Int e^(-W/z)`) is ARGUED [S; TRM s.2.3].
5. Unit and dilaton leaf: (2-32) equals `M_0(z/kappa)^-1 1^` [V through `z^8`; for all orders it follows from step 4 and the frozen unit].
6. Normalisation: `omega_{g,n} = Lambda^(g-1) Sum kappa^-(3g-3+n-Sum d) <..>^frozen prod d xi`; for `kappa = 1`: `(-5)^(g-1)`, no psi rescaling
   [W; TRM s.2.4]. Pairing with legs and D5 gives `<Phi_1..Phi_n>_h = (-5)^(h-1) <q_Phi1 .. q_Phin>_h` [W; LL s.1].
7. Free energies: dilaton equation on the CohFT side (`<tau_1(1)>_h = (2h-2) F_h`) and on the TR side (`(2-2h) F_h^TR = <Phi>_h`,
   `dPhi = y dx`), with `q_Phi = -tau_1(1)` (LD at `(0,1)`) [W; SYN2 M1 and the LD box].

## `lem:trm-xfun` (M6)
1. The linear loop equation [Q; [DW] EO Thm 4.4, (4-23)-(4-24); the sources cite only "s.4"].
2. Passage from (1) to "`Sum_c Res_c f(x) omega = 0` in one slot" [M; no source writes it] (suggested route in s.7 item G3).

## `prop:trm-legdict` (LD) [all W; SYN2 s.2 step 7.2 box]
1. Set-up: `D = d/dx = p^2/(p^5 - v) d/dp`; `xif_a`, `xif_{a,d} = D^d xif_a`; legs `Phi_0^b`; `D Phi_(k+1) = -Phi_k`.
2. Branch independence: each leg is analytic near every `c` for any branch of `log p`; a change of branch adds `2 pi i f(x)`; repeated
   integration by parts reduces `Res_c f(x) d xif_{a,d}` to `+-Res_c f^(d+1)(x) xif_a dx`, zero because `xif_a dx = eps_a v^(alpha_a) p^(m_a - 2) dp`
   is holomorphic at `c`.
3. (i) For `d >= 1`, integration by parts at each `c`: `<Phi, d xif_{a,d}> = -<D Phi, d xif_{a,d-1}>`; hence
   `<Phi_(k+1), d xif_{a,d}> = <Phi_k, d xif_{a,d-1}>`.
4. (ii) For `k >= 1`: `<Phi_k, d xif_a> = Sum_c Res_c Phi_(k-1) xif_a dx = 0`.
5. (iii) For `d >= 1`: `<Phi_0^b, d xif_{a,d}> = -Sum_c Res_c D^d(xif_a) dPhi_0^b`; this form is meromorphic on P^1, holomorphic at `p = 0`
   (`D` maps functions holomorphic at 0 to `O(p^2)`, `dPhi_0^b = O(p^-2) dp`) and `O(p^(1-4d)) dp` at oo; so the residues at the `c` sum to 0.
6. (iv) `<Phi_0^b, d xif_a> = -eps_a kappa_b v^(..) (1/5) Sum_c c^(m_a + n_b - 4)`, `dPhi_0^b = kappa_b v^(..) p^(n_b) dp`; non-zero only for
   `m_a + n_b = 4 mod 5`, i.e. `a = b`; there it equals `(1, -1, 1, 3, 2) = -c_b`, independent of `v` by homogeneity [W-: the table
   `(kappa_b, n_b) = (1,0), (1,-1), (-v,-2), (3,2), (2,1)` and `Sum_c c^N = 5 v^(N/5) [5|N]` are left to the reader; [DW] re-derived the five
   values by hand: they agree].
7. Combining (i)-(iv): `<Phi_k^b, d xif_{a,d}> = -c_b delta_ab delta_dk`; then D5 gives `q_{Phi_k^b} = -c_b tau_k(H^b)` [W].

## `lem:trm-lemG`, `lem:trm-resid`, `prop:trm-dii` (D(ii)) [all W; DR s.3]
1. `B_n := Sum_{j=0..n} x^j Phi_(n-j)/j!` (`= [z^n] e^(zx) Sum_k z^k Phi_k`).
2. `dPhi_m/dp = -x' Phi_(m-1)` (`m >= 1`) makes the sum telescope: `dB_n/dp = x^n Phi_0'/n!`, a Laurent polynomial.
3. `B_n in M` (`M` a `Q[v,x]`-module, `lem:trm-closure`); write `B_n = A + log p . B`, `B in Q[v,x]`.
4. `dB_n/dp = A' + B/p + log p . B'` has no log term, so `B' = 0`: `B` is a `p`-constant, equal to the `p^-1` coefficient of
   `x^n Phi_0'/n!` (`A'` has no `p^-1` term), i.e. `B = R_n := Res_{p=0} W^n dPhi_0/n!`.
5. Inverting, `Phi_k = Sum_j (-x)^j/j! B_(k-j)`; taking log coefficients gives Lemma G.
6. Residue identity: GI Lemma V gives `I^C_m = Res_{p=0} p^m e^(-W/z) dp`, so `[z^-n] I^C_m = (-1)^n Res p^m W^n dp/n!`; with
   `Y^C = (I_0, -I_-1, -I_3, I_2, I_1)` and `PiC(tau_n H^b) = (-1)^(n+1) [z^-n] Y^C_b`: `PiC(tau_n H^b) = -s_b Res p^(m_b) W^n dp/n!`,
   `(m_b, s_b) = (0,1), (-1,-1), (3,-1), (2,1), (1,1)`.
7. Compare with `dPhi_0^b = dp, dp/p, -v p^-2 dp, 3p^2 dp, 2p dp`; for `b = 2` use `Res W^n (p^3 - v p^-2) dp = Res d(W^(n+1))/(n+1) = 0`.
   Result `R^b_n = c_b PiC(tau_n(H^b))`, `c = (-1, 1, -1, -3, -2)`, every `n` ([DW] re-derived the five constants).
8. With `Pi_i(tau_k H^b) = PiC(tau_(k-i) H^b)` (MSG s.1) this is D(ii) [W; DR s.3]; "log-free iff `s in K_0`" because the `x^i` are
   independent [W; DP s.1].
* Inputs: GI Lemma V (PROVED there, GI s.2); the Stage-1 identification `Y^C_b = -(1/5) eta(S(z)H^2, H^b)` (MSG s.1, by code; signs
  cross-checked S01 C4) -- needed only if `K` is defined through `S(z)` (s.7 item G5).

## DP Lemmas 0-3 and Theorem D(iii) [all W; DP s.2-s.5, DR s.2 line-by-line]
* `lem:trm-count`: residues of the weights `{1,2,3,4}` of `b = 0,4,3,2` mod 4; the class-0 candidate `b = 2` needs `w - 5m >= 4` [W; DP s.2].
* `lem:trm-classL`:
  1. (a) `q = x - p^4/4`, so `Q[p,q] = Q[p,x]`, `v = pq = px - p^5/4` [W].
  2. `Q[v,x][T]/(T^5 - 4xT + 4v) -> Q[p,x]` is onto; its source is `Q[x,T]` (eliminate `v = (4xT - T^5)/4`), so it is an isomorphism;
     this gives the basis `1..p^4` [W-; DR Correction 6: add that `p, x` are algebraically independent].
  3. Converse of (a): for homogeneous `F in Q[p,q]` with most negative term `v^m p^-n` (`n > 0`, `m >= n`), `F - c v^(m-n) x^n in Q[p,q]`
     has a strictly higher most-negative exponent (lowest term of `x^n` is `v^n p^-n`); induction [W].
  4. (b) top-coefficient induction on the `x`-degree (works for `a_n in Q(v)`) [W]; (c) from (a), (b) [W]; (d) onto, kernel `Q[v,x]` by (b);
     `Q[v][p] cap Q[v][x] = Q[v]`; dimension count [W].
* `lem:trm-closure`: (i) compute `-Int v^c p^(a-c)(p^3 - v p^-2) dp` [W]; (ii) `d(log p . G(x)) = G dp/p + log p . g dx` and (i)-type
  integration [W; DR: both sides equal `log p (G - g_0) - Sum_(j != 0) g_j p^j/j`]; (iii) from (i), (ii); for `g in Q[v,x]`,
  `g(x) x' = d/dp G(x(p))` has no `p^-1` term [W].
* `lem:trm-modv`:
  1. Source basis of `F/vF` in weight `w`: one element if `4 does not divide w`; `{tau_(n-1)(H^2), tau_n(H)}` if `w = 4n >= 4`; `{tau_0(H)}` if `w = 0` [W].
  2. Target `N/vN = Q[p,q]/(Q[v,x] + pq Q[p,q]) (+) log p . Q[x]`; `Q[p,q]/(pq)` has basis `1, p^a, q^c`; the image of `Q[v,x]` is spanned by
     `1` and `x^n = p^(4n)/4^n + q^n` mod `pq`; so the basis is `{p^w}`, `{p^(4n), log p . x^n}`, or `{log p}` [W].
  3. `I` descends to `N/vN` (Lemma 2(iii)) [W].
  4. Images mod `v`: `Phi_k = (-1)^k p^(j+4k)/prod_(i=1..k)(j+4i)` (`j = 1,2,3`); `Phi_k^2 = (-1)^(k+1) p^(4k+4)/(4 prod_(i=1..k)(4+4i))`
     (using `q = -p^4/4` mod `Q[v,x]`); the `log p` component of `Phi_k^1` mod `v` is `(-1)^k x^k/k!` [W].
  5. In weight `4n` the matrix is triangular with non-zero diagonal; otherwise `1 x 1` non-zero [W].
* `thm:trm-lattice`:
  1. (a) = `lem:trm-count` [W].
  2. (b): `Lambda(s) in M` (Lemma 2), log-free by D(ii) (`s in K_0`), so `Lambda(s) in Q[p,q] in L` (Lemma 1(a)) [W; DP s.4].
  3. (d): `Lambda` graded `Q[v]`-linear, pieces finite-dimensional, zero in negative weight; `C = coker Lambda`; `C/vC = coker(Lambda mod v) = 0`
     (Lemma 3, right exactness); graded Nakayama (minimal `w` with `C_w != 0` gives `C_w = (C/vC)_w = 0`) [W]; `dim S_w = d(w) + e(w) =
     dim N_w` (Lemma 0, Lemma 1(d)), so bijective; `N` free [W; DP s.5; DR: `v` is a non-zero-divisor on `N`].
  4. (c): image in `Q[p,q]` by (b); injective (`Lambda(s) in Q(v)[x]` forces `Lambda(s) in Q[v,x]` by Lemma 1(b), zero in `N`, `s = 0` by (d));
     onto (`s := Lambda^-1(y (+) 0)` is log-free, so in `K_0` by D(ii)); dimensions `d(w)` [W; DP s.5, with O2/DR Corr. 2 for the
     dimension of `K_0 cap S_w`].
* `cor:trm-split`: by Lemma 2(ii) the log component of `Phi_k^1` is `(-1)^k x^k/k!` + lower `x`-degree; the log map on the `b = 1` triples is
  triangular onto `Q[v,x]_w`, hence bijective; its kernel is `K_0` (D(ii)) [W; DP s.5].
* `cor:trm-Ksplit`: by D(ii) `s in K` iff the log component of `Lambda(s)` is `x`-free, i.e. equals `Pi_0(s) in Q v^(w/5)`; write
  `s = s_0 + t` (Cor. 4); the leg of `v^(w/5) tau_0(H)` is `v^(w/5) log p`; injectivity on the `b = 1` span gives `t = Pi_0(s) tau_0(H)` [W; DP s.9].
* `prop:trm-coeff`: RD Lemma 0 (frames, by elimination from the top power; RD s.1.3) [W; RD]; the exponent of `f_j` from homogeneity
  (`w_a = c_a + 1`; `c_a + 1 >= 1` for `a != 1`, so the integer exponent is `>= 1`) [W-; SYN2 M5, one clause]; `kappa(X,k)` with `Q[v]`
  coefficients (MSG s.1(i)) and the count `d(w)` (D(iii)) give the basis [W; DR Correction 5(iii)].

## `thm:mlk-lreg` (L-regularity) [XY s.2.1, s.3; XYR s.3]
1. ABDKS Thm 7.1 gives `omega_{g,n}` by the closed formula (7.9) for `S = CP^1`, `y = z`, rational `x` with simple critical points,
   initial form `-y dx` [Q].
2. Each term of (7.9) is a Laurent polynomial in the `z_i` over powers of `q(z_i) = z_i^2 x'(z_i)` and `(z_i - z_j)`; the total has no diagonal
   poles and is regular at `z_i = 0` (`V != 0`); hence `omega_{g,n} = N(z; params) prod dz_i/prod q(z_i)^K` with `N` polynomial (on the line:
   `q = 4z^5 - v`, `K = 6g - 4 + 2n`) [W; XY s.2.1].
3. For legs in `L` the pairing is `(-1)^n Res_oo .. Res_oo prod Phi_j omega_{g,n}` (no poles at 0; functions of `x` drop out: `lem:trm-xfun`) [W-; XY s.3].
4. `1/q(z)` expands at `oo` with coefficients polynomial in the parameters (`z_i` and `z_i - z_j` are monic), so the pairing is a polynomial [W; XY s.3; XYR s.3].
5. Kernel sign: `omega^ABDKS_{g,n} = (-1)^n omega^engine_{g,n}` [W; XYR s.2, "PROVED here" + VERIFIED]. Transfer to the form D1 [M; s.7 item W5].

## `thm:mlk-gap` (Gap Theorem)
Route 1 [TRM s.4.1, all W after the quote]:
1. Hock Cor. 1.2, (1.3): `F_g = (1/(2-2g)) Sum_b Res_{z=b} [hbar^(2g-1)] dz Sum_{m>=2} (d^(m-1) z/dx^(m-1)) [u^m] exp[(S(u hbar d_z) - 1) x u]/(u hbar)`,
   `S(t) = sinh(t/2)/(t/2)`, `b` over the points where `y dx` or `x dy` has a non-zero residue (here `b = a`, `b = oo`) [Q].
2. Extraction: `F_g = (1/(2-2g)) Sum_b Res I_g dz`, `I_g = Sum_{j=1..g} (1/j!) D^(2g-2+j)(z) Sum_{k_1+..+k_j = g} prod_i x^(2k_i)/(4^(k_i) (2k_i+1)!)` [W].
3. At `z = oo`: every term has degree `<= 1 - (2g-2) r - 2g <= -5`, residue 0 [W].
4. At `z = a` (`= 0`): `D^n(z) = (-1)^n n! z^(n+1) v^-n (1 + O(z^2))`, `x^(2k) = (2k)! v z^-(2k+1)(1 + O(z^(2k+1)))`; lowest order exactly `z^-1`, the
   regular part of `x` drops out: `Res = v^(2-2g) A_g`, `A_g = Sum_j (-1)^j (2g-2+j)!/j! [t^g] phi(t)^j`, `phi = Sum_k (t/4)^k/(2k+1)` [W].
5. `A_g = (2g-2)! [t^g](1 + phi)^(1-2g)`; `1 + phi = artanh(s)/s`, `t = 4s^2`; Lagrange inversion gives `A_g = -B_2g/(2g)`, so
   `F_g = B_2g/(2g(2g-2)) v^(2-2g)` [W]. Sign remark: Hock's s.3.2 example prints `B_2g/(2g(2-2g))` (slip) [W; TRM s.4.1].
Route 2 [XY s.2.3; XYR s.2]: monomial count near 0 and oo (each `T_r` regular at 0, `O(z^(5-10g))` at oo) [W]; Borel summation of
`Sum_r (r-1)! [t^(r+1)] e^(E_0)` with `eps s = tanh(eps u)` and `(x/sinh x)^2 = Sum (1-2g) 2^(2g) B_2g x^(2g)/(2g)!` [W]; the polynomial
part of `Phi` contributes 0 (negative weighted degree) and the `-V log z` part gives `Phi_g(0) - Phi_g(oo)` (contour lemma; XYR two-line
proof: `Res_a log z dPhi = -Res_a Phi dz/z`) [W]. Frozen form: multiply by `(-5)^(1-h)` (`thm:trm-ident` step 7) [W; SYN2 M8].
* `lem:mlk-vanper`: route 2 at `g >= 1` [Q + W; XY s.2.3]; orientation `<omega_{g,1}, -log z> = -Int_oo^0 omega_{g,1}` [W; XYR s.2].
