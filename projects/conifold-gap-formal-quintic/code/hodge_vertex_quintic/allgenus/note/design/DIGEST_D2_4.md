# DIGEST_D2, part 4: s.2C (QRR) statements; s.3 Proofs, part A

Continuation of `DIGEST_D2_3.md`. Proofs of part B (formal lemmas) and s.4 conventions: `DIGEST_D2_5.md`.

## 2C. The input (QRR)

**QRR-1** `thm:qrr-input` (thm, quoted input; proposed placement: the statement of inputs in s.1 and the set-up). *Statement.* In
Givental's convention (`S ~ R e^(u/z) C`, `R = 1 + O(z)`; legs `R^-1(psi)`, edges `(eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T)/(psi + psi')`,
LIT s.1.2), the R-matrix of the formal-quintic CohFT satisfies `R(q = 0) = diag(h_(w_0), .., h_(w_4))`, to all orders in `z` (its `z^(5k)`
part included). *Form used:* ONLY `A^QRR_n(q = 0) = h_n` (H9), i.e. the flat colour-0 QRR column has `Phi^QRR(z, u = 1) = h(z)`
("the u = 1 convention of ref_rows" of EVID, folded into (QRR) by RF-C4) (M5 s.2; SYN5 s.1; RS s.3 row (QRR), PS s.4.3).
*Sources* (LIT s.1.2; SYN5 s.1; M5 s.2):
* Guo-Janda-Ruan arXiv:1812.11908v1, eq. (17), in the proof of Lemma 5.2: `C_alpha(z) = exp(Sum_(k>0, beta != alpha) b_k z^(2k-1)/
  (zeta^alpha lambda - zeta^beta lambda)^(2k-1)) exp(Sum_(k>0) b_k z^(2k-1)/(-5 h_alpha)^(2k-1))`, `h_alpha = zeta^alpha lambda`, with the Birkhoff
  factorisation `S_(i alpha)(z) C_alpha(z)^-1 = e^(u_alpha/z) R_(i alpha)(z)`; hence `R(q = 0)_(alpha alpha) = C_alpha^-1 = h_(w_alpha)`. GJR's
  "lambda-twisted theory" (s.3.1, eq. (9): `t = 0`, `lambda_i = zeta^i lambda`) is the formal quintic; s.4.1: a semisimple CohFT
  reconstructed by Givental-Teleman. Read in the worker's own PDF text extraction (LIT3 s.5).
* Guo-Janda-Ruan arXiv:1709.07392v1, eqs. (19)-(20), in the proof of Lemma 4.6, at `t = 0` (twist weight `t - 5 lambda_i`; Lemma 4.6
  writes the `q = 0` diagonal of `R` through `z^3`; base weights `lambda_i = xi^i lambda`).
* Both justify it only by citing Givental arXiv:math/0008067 Thm 1 with the Remark after it, and Coates-Givental arXiv:math/0110142
  Thm 1, Cor. 1 (LIT s.1.2-s.1.3). Neither paper proves it (LIT s.1.2 "Status").
*Status.* Quoted (M5 s.2; SYN5 s.1). *Conflict settled* (SYN5b K9, O11; RF2 s.5; CP4 s.5 item 2): M4 s.2 cited arXiv:math/0108100 Thm 9.1
(the untwisted equivariant statement, cotangent characters) and arXiv:1702.06096 eq. (5) (the definition of the theory, no
`R(q = 0)`); the correct citation is GJR. LP's `R_k` (arXiv:1803.01409 s.2.4.4 eq. (18), Prop. 5) carry the Zagier-Zinger
normalisation and vanish at `q = 0` for `k >= 1`: not the QRR-normalised R-matrix (LIT s.1.1).

**QRR-2** `rmk:qrr-sign` (rmk). *Statement.* The twist weight is `-5 zeta^i`: in `h_(w_i)` the O(5) factor is
`exp(+Sum_k b_k (z/(5 lambda_i))^(2k-1)) = exp(-Sum_k b_k z^(2k-1) (-5 lambda_i)^(1-2k))`, the tangent factor is
`exp(-Sum_k b_k Sum_(j != i) (z/(lambda_i - lambda_j))^(2k-1))`. Settled three ways (LIT s.1.4): (a) LP arXiv:1803.01409v2 eq. (26) with
Mumford's relation and Givental arXiv:math/0108100v2 Thm 4.1 (QRR-3); (b) Coates-Givental arXiv:math/0110142v2 Cor. 1 (an Euler-class
twist with Chern root `lambda + rho` contributes `exp(+Sum_m b_m (z/(lambda + rho))^(2m-1))`; at `p_i` the root of O(5) is `5 lambda_i`); (c)
GJR (17). Pinned independently by `h_1 = ph_1 = -3/20` (`ph_1` is forced by HPF; the weight `+5` gives `-11/60`) (RF2 s.5; CP3 s.3.6).
*Status.* "PROVED given the quoted statements, whose signs need a human reading" (SYN5b O15: GJR and CG were read in a text extraction
that drops minus signs; signs reconstructed from ar5iv or the structure of the formula). VERIFIED: from LP (26) alone,
`N_(2,0) = -19/9600`, `N_(3,0) = 1/224000` (the project's localisation values); weight `+5` gives `-97/28800`, `1/81000` (LIT s.1.5).
SYN5b K14: no conflict between CP3 s.3.6, LIT s.1.4, Stage-4 (H) at `n = 1`, GEN T2-T3.

**QRR-3** `rmk:qrr-vertex` (rmk), the vertex-level derivation (LIT s.1.4 (a); PROVED given its quotes, M5 s.2). LP (26):
`H_g^(p_i) = [e(E^* (x) T_(p_i) P^4)/e(T_(p_i) P^4)] [(5 lambda_i)/e(E^* (x) (5 lambda_i))]`. With `Lv(w) := e(E^* (x) L_w) = Sum_k (-1)^k
lambda_k w^(g-k)` and Mumford's `c(E) c(E^*) = 1` in the form `Lv(w) Lv(-w) = (-1)^g w^(2g)`: `H_g^(p_i) = W_i^(g-1) prod_(a in w_i) Lv(a)/a^g`,
`W_i = prod_(j != i) (lambda_i - lambda_j)/(5 lambda_i)` (the genus weight, not part of R). Then `Lv(a)/a^g = exp(-Sum_k (2k-2)! ch_(2k-1)(E)
a^(1-2k))` (`ch_2k(E) = 0`, `k >= 1`), and Givental's Thm 4.1 with `s_k = -(2k-2)! Sum_(a in w_i) a^(1-2k)` gives the operator of `h_(w_i)`.
Sources quoted there: Mumford's relation (FP arXiv:math/9810173 s.2.1, unnumbered; Mumford 1983 only as quoted); FP Prop. 1-2
(Prop. 2 = operator form of Thm 4.1). *Gap stated by LIT s.1.4:* that the R-matrix of the full CohFT restricts at `q = 0` to these
vertex R-matrices is quoted, not proved: the degree-0 part is the direct sum of the five vertex theories
(`Mbar_(g,n)(P^4, 0) = Mbar_(g,n) x P^4`) and the R-matrix of a semisimple CohFT is unique (QRR-4).

**QRR-4** (uniqueness; for a `\draftnote`). Teleman arXiv:0712.0160v3: Theorem (2) in s.1.9 (CohFTs with a given semisimple Frobenius
algebra form the orbit of the trivial theory) and Prop. 8.5 (ii) in s.8.2 (inhomogeneous ambiguity: theories over a given Frobenius
manifold differ by a Hodge-type twist), read ONLY through the ar5iv rendering and a summariser (SYN5 s.4.3). "The step 'the R-matrix
is fixed by the full CohFT, so `R(q = 0)` is the vertex R-matrix' is not literally stated in what I saw: a human check" (SYN5 s.4.3;
SYN5b s.7 item 5: does the full CohFT fix the Hodge twist of Prop. 8.5 (ii)?).

**QRR-5** (what (QRR) does and does not pin; for s.1 of the note). Theorem Gamma computes `c(z; 1) = h(z)` with no reference to the
actual R-matrix, so "`D = 1`" is proved for the R-matrix normalised by (QRR) (RF2 s.5). The overall sign and the weight `-5` are pinned
by `h_1 = ph_1`; the `z^(5k)` parts are not pinned by any check; the genus-3 code test `F_3(0) = 1/224000` from the project's graph sum at
`u = 1` is optional and NOT done (SYN5b s.7; RS3 s.8; LIT s.1.5: `N_(3,0)` there is a localisation number, not that test).

## 3A. Proofs of part A (each step: where it is written out)

**Proof of GAM-1** (all written out). (1) `S_0`, `A_St` holomorphic on `Delta`, `S_0'(L) = 0`, `S_0''(L) = 5/(L(1 - L^5)) != 0` (CP s.1.6).
(2) By (F1), `c_n(L)` is a polynomial over Q in the Taylor coefficients of `S_0` (orders `>= 3`), of `S_1 - S_1(L)` and of the `S_2k`, and
in `-1/S_0''(L)` (CP s.1.6). (3) These are the explicit expressions of D6, polynomials in `L^-1` and the `P_k in (1 - L^5)^-k Q[L]` (CP s.1.6;
RF s.1 checked the generating function by hand). (4) Hence `c_n in Q[L, L^-1, (1 - L^5)^-1]`. (5) Choices: Taylor coefficients of order
`>= 1` carry no branch; `e^(S_1(sigma) - S_1(L))` is the branch equal to 1 at `L`; `lambda` is forced by criticality; other branches shift
the phase by constants (CP4 s.4 (B1), (B2), (B5), (B7)); RA s.1.3 last item. Complete.

**Proof of GAM-2.** Route RF-C5 (RF s.1, s.1.7 item "A shorter route"): (1) the second derivatives of the two flat phases coincide as
functions; (2) the normalised amplitudes coincide (odd exponents); (3) `L` is critical in both by the choice of the logarithm; (4) `Nhat`
depends only on these data (CP s.1.7; RA s.1.3). Written in one sentence each. *Missing sentence:* the cancellation of the `(log z)/z`
terms for POLY's form (coefficient `-5 sigma + Sum_j (sigma - zeta^j) = 0`) is not written for that form; it follows from the reflection
route (CP s.1.7), which reduces POLY's integrand to GAMMA's times `e^(c/z)` and a linear phase. Route CP s.1.7: reflection identity
(SP-REFL, RL s.3.3), `log y - log(-y) = i pi eps_j`, `e` even, `e = 0` (written in CP s.1.7; checked RF s.1, RA s.1.1).

**Proof of GAM-3** (RA s.3.2; re-derived LF s.2, LA s.1.4; all written out).
1. `d_s = P'(s) + alpha0 m`, `d_m = alpha0 s` vanish at 0; `H = [[2 p_2, alpha0], [alpha0, 0]]`, `det H = -alpha0^2`;
   `C = -H^-1 = [[0, -1/alpha0], [-1/alpha0, 2 p_2/alpha0^2]]`: `C_ss = 0`, `C_sm = -1/alpha0`, `C_mm = 2 p_2/alpha0^2` (LF s.2).
2. `Phi_3 = Sum_(k>=3) p_k s^k` (the term `alpha0 s m` is quadratic). Scaled variables `s = z^(1/2) eta_s`, `m = z^(1/2) eta_m`: the integrand is
   `e^(J eta_m) Pi(eta_s)`, `J = c z^(1/2)`, `Pi = exp(Sum_(k>=3) p_k z^(k/2 - 1) eta_s^k) A(z^(1/2) eta_s, z)`, a function of `eta_s` alone.
3. Moments: since `C_ss = 0`, each `eta_s` pairs with a distinct `eta_m`: `<eta_m^r eta_s^a> = r!/(r-a)! (r-a-1)!! C_sm^a C_mm^((r-a)/2)` for
   `r - a` even `>= 0`, else 0 (LF s.2; "`= 0` for `r < a`" implicit in RA). With `(2j-1)!!/(2j)! = 1/(2^j j!)`:
   `Sum_r J^r/r! <eta_m^r eta_s^a> = (J C_sm)^a e^(J^2 C_mm/2)`.
4. Hence `<e^(J eta_m) Pi> = e^(J^2 C_mm/2) Pi(J C_sm)`, termwise in `z^(1/2)` (finitely many terms per order).
5. `z^(1/2) J C_sm = -c z/alpha0 = delta`, `J^2 C_mm/2 = p_2 delta^2/z`, `Pi(J C_sm) = exp(Sum_(k>=3) p_k delta^k/z) A(delta, z)`: `N = exp(P(delta)/z) A(delta, z)`.

**Proof of GAM-4** (RA s.3.3, with LF s.1, s.3 and LA s.1.1-s.1.4; all written out once LF-C1 is inserted).
0. *Definition* (LF s.1, LA s.1.1). `S'(sigma) = lambda + g'(sigma) + alpha0 log a0(sigma)`, so `lambda(L)` makes `L` critical;
   `S''(L) = g''(L) + alpha0^2/a0(L) = (g''(L) a0(L) + alpha0^2)/a0(L)`, non-zero on `U*` (LF-C3). The Taylor coefficients of order `k >= 2`
   of `a0 log a0`, `(-1)^k alpha0/(k(k-1)(L - s*)^(k-1))`, carry no branch: (i). By (F1) every `c^F_n` is holomorphic on `U*`, a priori with
   a pole at `s*` (LA s.1.1).
1. *Two variables.* `Omega_2 := exp([lambda sigma + g(sigma) + a0(sigma) nu - e^nu]/z) e^(M0 nu) F(sigma, z) d sigma d nu`, a FORMAL integrand
   (LA-C1: the actual Euler integral of `Gamma(M0 + a0/z)` carries `z^(M0 + a0/z)`, whose `sigma`-dependent part would shift `lambda` by
   `alpha0 log z`; `c^F` is defined without it). At fixed `sigma` near `L`: `d_nu = a0 - e^nu`, critical point `nu_c(sigma) = log a0(sigma)` (the
   branch used in `lambda`, continuous in `sigma`), Hessian `-a0(sigma) != 0`. By (L1a) with `chi = exp` (`J = e^nu`, `J(nu_c) = a0`) and (L3)
   (`a = a0`, `M = M0`, measure `dx`): `N[a0 nu - e^nu, e^(M0 nu); nu_c] = a0 N[a0 log x - x, x^(M0 - 1); a0] = a0^M0 Poch_M0(z/a0) exp(Sum_k b_k
   (z/a0)^(2k-1))`; `F` factors out (`C[[z]]`-linearity, LA-C2). Outer phase `S` ((L1b)(i)); total Hessian `[[g''(L), alpha0], [alpha0, -a0(L)]]`
   (`a0'' = 0`), `det = -g''(L) a0(L) - alpha0^2 = det H_vv . S''(L) != 0`; `beta = (a0(sigma)/a0(L))^(-1/2)`, `beta B = a0(L)^(1/2) A` (the
   constant is harmless, LF s.3.1). By (L1b): `Nhat[Omega_2; (L, nu_c(L))] = c^F(z; L)` for `L in U*`.
2. *Chart.* Translate `nu = mu - lambda/alpha0` (Jacobian 1). Because `a0` is EXACTLY `alpha0 (sigma - s*)`:
   `lambda sigma + a0 nu = lambda s* + alpha0 (sigma - s*) mu`; `e^nu = Q e^mu`, `Q := e^(-lambda/alpha0)`; `e^(M0 nu) = Q^M0 e^(M0 mu)`. Phase
   `lambda s* + Psi_Q`, `Psi_Q := g(sigma) + alpha0 (sigma - s*) mu - Q e^mu`; amplitude `Q^M0 e^(M0 mu) F`. (F2) removes `lambda s*`; `Q^M0`
   drops out of `Nhat = N/N_0` (LF-C2). *Q is a coordinate:* `Q(L) = a0(L) e^(g'(L)/alpha0)` is single-valued (a change of branch of `log a0`
   changes `lambda` by `-2 pi i k alpha0`), holomorphic at `s*`, `Q(s*) = 0`, `Q'(s*) = alpha0 e^(g'(s*)/alpha0) != 0` (LF s.3.2). *The family:* at
   `Q = 0` the critical point is `(s*, mu_*)`, Hessian `[[g''(s*), alpha0], [alpha0, 0]]`, `det = -alpha0^2`; by the holomorphic IFT (CP s.1.4)
   there is a holomorphic non-degenerate family `p(Q) = (sigma_c(Q), mu_c(Q))`, `|Q| < r`; the critical equations `g'(sigma) + alpha0 mu = 0`,
   `alpha0 (sigma - s*) = Q e^mu` give `sigma_c(Q) = s* + Q e^(mu_*)/alpha0 + O(Q^2)`, `sigma_c'(0) != 0`. *Identification (LF-C1, LF s.3.2):* for
   `0 < |Q| < r` put `L := sigma_c(Q)` (`!= s*`, since `alpha0 (L - s*) = Q e^(mu_c) != 0`), choose a branch of `log Q`, put
   `lambda := -alpha0 log Q`, `nu_c := mu_c + log Q`; then `e^(nu_c) = a0(L)` and the `sigma`-equation gives `lambda = -g'(L) - alpha0 nu_c`, the
   `lambda` of step 1 for this branch; the translation maps `(L, nu_c)` to `p(Q)`. So `Nhat[Psi_Q, e^(M0 mu) F; p(Q)] = c^F(z; sigma_c(Q))`, and
   since `sigma_c` maps the disc biholomorphically onto a neighbourhood of `s*`, `c^F_n = chart_n o sigma_c^-1` there. No uniqueness
   statement is needed (LF s.3.2; LA s.1.3 gives the converse identification). *Holomorphy:* `Psi_Q` and `e^(M0 mu) F` are holomorphic in
   `(sigma, mu, Q)` near `(s*, mu_*, 0)`; by (L1c)(i) every `N_n(Q)` is holomorphic; `N_0(Q) = e^(M0 mu_c(Q)) F(sigma_c(Q), 0) != 0` near 0 (uses
   `F(s*, 0) != 0`), so `Nhat_n(Q)` is holomorphic and at `Q = 0` equals `Nhat_n` of the `Q = 0` data: (ii). (LF s.3.2: explicitly the chart
   coefficients are polynomials in Taylor data, in `Q e^(mu_c) = a0(L)` and in `-H^-1`, with `1/det H = -1/(g''(L) a0(L) + alpha0^2)`.)
3. *Value.* With `s = sigma - s*`, `m = mu - mu_*`: `Psi_0 = g(s*) + gt(s) + alpha0 s m` exactly (the linear terms `g'(s*) s + alpha0 mu_* s`
   cancel; no pure `m`-term, no `s^k m` with `k >= 2`); `e^(M0 mu) = e^(M0 mu_*) e^(M0 m)`. GAM-3 with `P = gt`, `c = M0`, `A = F(s* + s, z)`:
   `N = e^(M0 mu_*) exp(gt(delta)/z) F(s* + delta, z)`; divide by `N_0 = e^(M0 mu_*) F(s*, 0)`: (iii). □
*Removable singularity* (LA s.2.1; FSL3 s.8 warning 2): elementary; for the quintic `c_n` is rational (GAM-1), so it simply has no pole.
*The `M0 != 0` question* (LF s.3.3): the Pochhammer factor enters only in step 1, as an exact identity; at `Q = 0` the same `M0` reappears
as the source `e^(M0 m)`, which GAM-3 turns into the shift `delta`.

**Proof of GAM-5.** (ii) LMI at `s* = 1` (LF s.4.1; LA s.1.5): `alpha0 = -1`, `a0 = 1 - sigma`, `M0 = 0`; the `j = 0` factor
`Gamma((1 - sigma)/z)` contributes phase `a0 log a0 - a0`, amplitude `a0^(-1/2) exp(Sum_k b_k (z/a0)^(2k-1))`;
`g = [5 sigma log 5 sigma - 5 sigma] + Sum_(j>=1) [(zeta^j - sigma) log(zeta^j - sigma) - (zeta^j - sigma)]`;
`F = (5 sigma)^(1/2) prod_(j>=1) (zeta^j - sigma)^(-1/2) exp(Sum_k b_k z^(2k-1) [(5 sigma)^(1-2k) + Sum_(j>=1) (zeta^j - sigma)^(1-2k)])` (LF-C4 =
LA-C4: `F` is `e^(S_1)` without its `j = 0` factor, so it contains `(5 sigma)^(1/2)`). Then `S = S_0`, `A = A_St`, LMI's `lambda` is CP's, so
`c^F = c` near 1; `g` is holomorphic at 1 (`5 != 0`, `zeta^j - 1 != 0`), `F(1, 0) != 0`; `Q = e^lambda = -q`. LMI (ii) + GAM-1: no pole at
1. Value (`delta = 0`): `F(1, z)/F(1, 0) = exp(Sum_k b_k T_k z^(2k-1)) = h(z)` (D3). (iii) CP3 s.3.9 (written out; RF2 s.4 checked every
step): `P_k(zeta L) = zeta^-k P_k(L)`; the Taylor data at `zeta L` are `zeta^-k`, `zeta^-k`, `zeta^(1-2k-r)` times those at `L`; `t = zeta t'`
((L1a), constant Jacobian); the phase becomes `zeta` times the phase at `L`, the normalised `e^(S_1)` part is unchanged, `z^(2k-1) S_2k`
becomes `(zeta^-1 z)^(2k-1) S_2k(L + t')`; Lemma H with `mu = zeta`: `c(z; zeta L) = c(zeta^-1 z; L)`. Consequences: regular at `zeta^i`
with value `h(zeta^-i z)`; exponents `= 4n mod 5`; `c_n in Q[L, L^-1]`. Alternatively LMI at `s* = zeta^i`, `a0 = zeta^i - sigma`,
`alpha0 = -1` (LF s.4.1). (iv) GAM-9, GAM-8.

**Proof of GAM-6** (CP3 s.3.9, RF2 s.4). `q(zeta^i L) = q(L)`; the colour-`i` point `zeta^i L` of the parameter `L` is the colour-0 point
of the parameter `zeta^i L` (logarithm chosen as in CP4 (B1)); GAM-1 (dependence on the critical point only) and GAM-5 (iii).

**Proof of GAM-7** (CP3 s.3.10 with RF-C2; RF2 s.5 item 2). (1) `sigma = L/X` is analytic near `X = 1` ((L1a)); its Jacobian `-L/X^2`
combines with the factor `X^2` of `e^(S_1(L/X))` into POLY's `(1 - u X^5)^(-1/2)`: part of the hatted amplitude (RF-C2). (2)
`S_0(L/X)/z = L^4 S_0(L/X)/w`; Lemma H with `mu = L^-4` identifies the normalised `z`-series of the `sigma`-form with the normalised `w`-series
of the pulled-back form at `w = L^4 z`. (3) The pulled-back form is POLY's hatted integrand up to a constant (GAM-2); its normalised series
is `Phi^SP(w; u)`. (4) Both sides are rational in the parameter ((L1c)(i)); hence `Sum_n c_n(L) z^n = Sum_n ph^SP_n(L^-5) L^(4n) z^n`. (5) At
`L = 1`: `c_n(1) = ph^SP_n(1)`; `ph^SP_n` is regular at `u = 1` by GAM-5 (ii) or Theorem P(a).

**Proof of GAM-8** (LF s.4.3). GAM-1 + no pole at `L^5 = 1` (GAM-5 (ii)-(iii)) + no pole at `L = 0` (GAM-9) give `c_n in Q[L]`; the Z/5
exponents give the form. With GAM-7 the coefficient of `L^(4n-5k)` is `[u^k] ph^SP_n`; `4n - 5k >= 0` gives the top-degree bound.

**Proof of GAM-9** (LF s.4.2 with LF-C5; LA s.1.5, s.2.3). (1) Near `sigma = 0` only `Gamma(1 + 5 sigma/z)` degenerates: `s* = 0`, `alpha0 = 5`,
`a0 = 5 sigma`, `M0 = 1` (`(5 sigma)^(1/2) = a0^(M0 - 1/2)`, `Poch_1 = 1`); `g`, `F` come from the five factors `j = 0..4`, regular at 0;
`lambda` is CP's; `Q = (-q)^(-1/5)`; `delta = -z/5`. (2) LMI: regular at `L = 0`. (3) Value as a product: with `a_j = zeta^j - sigma`,
`e^(g/z) F = (2 pi)^(-5/2) z^(-5/2) z^(-5 sigma/z) prod_j Gamma_St(a_j/z)` (`Sum_j zeta^j = 0`); at `sigma = delta`: `a_j(delta)/z = zeta^j/z + 1/5`,
`z^(-5 delta/z) = z`, `e^(g(delta)/z) = e^(g(0)/z) e^(-g'(0)/5) e^(gt(delta)/z)`; so the value is `prod_j Gamma_St(zeta^j/z + 1/5)/[(2 pi)^(5/2) z^(3/2)
e^(g(0)/z) e^(-g'(0)/5) F(0, 0)]`. (4) SP-5118 at `w = zeta^j/z`, `theta = 1/5` (`log w = log zeta^j - log z`): the `log z` terms give
`z^(-5(theta - 1/2)) = z^(3/2)`; the `1/z` terms give `e^(Sum_j zeta^j log zeta^j/z) = e^(g(0)/z)`; the constants give
`exp((theta - 1/2) Sum_j log zeta^j) = e^(-g'(0)/5) F(0, 0)` (`g'(0) = -Sum_j log zeta^j`, `F(0, 0) = prod_j zeta^(-j/2)`, same branches).
(5) Everything cancels: `c(z; 0) = exp(Sum_(n>=1) (-1)^(n+1) B_(n+1)(1/5) z^n Sum_j zeta^(-jn)/(n(n+1)))`; `Sum_j zeta^(-jn) = 5 [5 | n]`;
`n = 5l` gives `Ginf`. All written out in LF s.4.2.

**Proof of GAM-10** (RA s.3.4; LF s.4.3). `[u^(4m)] ph^SP_(5m)` is the coefficient of `L^0` in `c_(5m)` (GAM-7), i.e. `c_(5m)(0)`, which
is `[z^(5m)] Ginf` (GAM-9). With P(b) (`ph^SP = ph^can`): `f^can_m = [z^(5m)] Ginf`.

**Proof of GAM-11** (CP3 s.3.1-s.3.6, all written out; RF s.3.1-s.3.6, RA s.1.3 CONFIRMED). (a1) The six Euler substitutions give
`z^E`, `E = -1 - (1/z)[5 sigma + Sum_j (zeta^j - sigma)] = -1` identically in `sigma` (s.3.1). (a2) Critical point, `H_7`, `det H_7` (s.3.2; the
rows multiplied by `L(1 - L^5)`, `5L`, `zeta^j - L` give a polynomial matrix of determinant 5). (a3) Order (a): (L1b) with inner `(y, x_j)`,
(L0), (L3) with `a = 5 sigma`, `M = 1` and `a = zeta^j - sigma`, `M = 0`; `beta B = const z^-1 A_St`: `Nhat_7 = c` (s.3.3). (b) Order (b): (L2)
with `w = (y, x_1..x_4)`, `d_(x_0) G = -1/x_0`, `X_0(w)`, `sigma_c = 1 - X_0`, `beta N = z^-1 (1 - L)^-1 prod_(j>=1) x_j^-1`, critical value
`lambda + psi(w)`: `Nhat_7 = Nch` (s.3.4). (c) Both orders compute the same `N` ((Fub) is an identity of series); `Nch_n` rational and
regular at 1 since `1/det H_ch = -L^2 Lam1(L)^2`; two rational functions equal on `Lam` are equal (s.3.5). (d) At `L = 1` the chart is a
product (L0) of the Euler integrals of `Gamma(1 + 5/z)` (`a = 5`, `M = 1`, `dy`) and `Gamma((zeta^j - 1)/z)` (`a = zeta^j - 1`, `M = 0`);
(L3) gives `exp(Sum_k b_k T_k z^(2k-1)) = h` (s.3.6). Branches: CP4 s.4.

**Proof of GAM-12** (CP3 s.3.7; RF s.3.7). (1) `psi_q = psi_0 + q K_1`, `K_1 = y^5/(x_1 x_2 x_3 x_4)`; (L1c)(iii) expands the SP at the moving
point coefficientwise in `q`. (2) Each `q^d` term factorises (L0) into Euler integrals with `M = 5d + 1` (`y`) and `M = -d` (`x_j`); with
`Gamma_St(w + 1) = w Gamma_St(w)` the ratio to `d = 0` is `q^d I_d(z)` (powers of `z` cancel; signs `(-1)^(4d) = 1`). (3) The prefactors of the
moving point give `L`. (4) `D_q Dpsi = q (d_q psi_q)(p_ch) = K = L - 1` (envelope theorem). All written out.

**Proofs of the remarks GAM-13..GAM-16.** *GAM-13* (LA s.2.3 Corollary): RA's proof verbatim with SP-L3C in step 1; `Q^M0 :=
e^(-M0 lambda/alpha0)` is a constant; GAM-3 allows any `c in C`. Written in two lines; refereed by no second agent (RS A3).
*GAM-14:* at `sigma = 1` the extra factor `5 sigma` of `Gamma(1 + 5s) = (5s) Gamma(5s)` is `z`-free and `delta = 0` (LA s.2.4), equivalently
`Poch_0 = Poch_1 = 1` in (L3) (RF s.1); at `L = 0` with `M0 = 0`, LMI's value is `F(0, z)/F(0, 0) = exp(Sum_k b_k z^(2k-1) Sum_j zeta^(j(1-2k)))
= exp(5 Sum_(5 | 2k-1) b_k z^(2k-1))` (LA s.2.4 states the value; the one-line evaluation is mine, from LMI (iii) with `delta = 0`; `\draftnote`).
*GAM-15 (a)* (GEN s.3.2 Lemma 4, written out there): near `sigma = lambda_i` the critical equation reads `log(+-Q) + l log(l (sigma - lambda_i))
- log(lambda_i - sigma) + O(1) = 0`, i.e. `(l - 1) log|sigma - lambda_i| = -log|Q| + O(1)`; for `l >= 2` the right side tends to `+oo`, for `l = 1`
it says `log|Q| = O(1)`; both impossible as `Q -> 0`. (An analytic statement about actual critical points; remark only.)
*GAM-15 (b)* (LF s.1, two sentences): (L1a) with `sigma - s* = (L - s*) t` and Lemma H with `mu = a0(L)` give `c^F(z; L) = C(z/a0(L))`; if
`C_n != 0` for some `n >= 1`, `C_n a0(L)^-n` has a pole at `s*`; LMI (ii) excludes it (and LMI (iii) gives the value 1 directly).
*GAM-16:* each case is an instance of LMI (LA2 s.5); the `z`-powers of a general summand cancel up to a `sigma`-independent factor
(GEN s.2 Lemma 3, written out, R0); the summand/I-function link is GEN Lemma 2 (residue ratios computed factor by factor, R0).
