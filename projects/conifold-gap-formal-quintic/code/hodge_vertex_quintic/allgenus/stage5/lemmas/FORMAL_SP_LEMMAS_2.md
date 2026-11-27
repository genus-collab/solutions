# Stage 5 / lemmas, part 2: FL1, translation, parameters (L1c, FL4), changes of variables (L1a)

Continuation of `FORMAL_SP_LEMMAS.md` (s.0-s.2: definitions, universality, Gaussian toolkit). Same notation and labels.

## 3. Exact forms and translations

### 3.1 FL1: the formal SP kills exact forms (PROVED)
**Lemma FL1.** Let `F, c, H, C` be as in s.1.4 and let `h_1, .., h_n` be amplitudes (in `A[[t]]((z))`, or any `h_i` with
`e^(F_3/z) h_i in R_n(A)`). Then
`SP^norm_c[ F ; Sum_i (d_i h_i + z^-1 (d_i F) h_i) ] = 0`,  i.e.  `SP_c[ d(e^(F/z) iota) ] = 0`  for `iota = Sum_i (-1)^(i-1) h_i dx_1..^dx_i..dx_n`.
The same holds for every sum of such expressions that converges in `R_n(A)` after multiplication by `e^(F_3/z)`.
*Proof.* Put `g_i := e^(F_3/z) h_i in R_n(A)`. `d_i` is a continuous derivation of `R_n(A)` with `d_i e^(F_3/z) = z^-1 (d_i F_3) e^(F_3/z)`
(termwise), so `d_i g_i = e^(F_3/z) (d_i h_i + z^-1 (d_i F_3) h_i)`. Since `d_i F(c + t) = (H t)_i + d_i F_3(t)`,
`e^(F_3/z) Sum_i (d_i h_i + z^-1 (d_i F) h_i) = Sum_i (d_i g_i + z^-1 (H t)_i g_i)`, whose `E_zC` vanishes by (G3). QED
*Remark (POLY s.6, Theorem P(a)).* There `g(sigma) - g(sigma + z) = -d/dsigma Sum_(k>=1) z^k g^(k-1)/k!` with `g = e^(Phi/z) A`. The sum
converges in the sense of FL1: `e^(-Phi/z) g^(k-1)` is a polynomial in the `Phi^(j)/z` and the derivatives of `A`, with
`sum of the j's <= k - 1`; since `Phi' = O(t)` has weight `>= 1/2` and `Phi^(j)` (`j >= 2`) weight `>= 0`, `z^k e^(-Phi/z) g^(k-1)` has weight
`>= k - (k-1)/2 = (k+1)/2`. So FL1 gives `SP[g(sigma) - g(sigma + z)] = 0`; Lemma T(b) below gives the same directly.

### 3.2 The translation lemma (PROVED)
**Formal implicit function theorem (IFT; PROVED).** Let `B` be complete and separated for an ideal `I` (`B = A[[y]]`, `I = (y)`; or
`B = A_0[[q]]`, `I = (q)`). Let `F in B[[t]]` with `grad F(0) in I^n` and `Hess F(0) in GL_n(B)`. Then there is a unique `c in I^n` with
`grad F(c) = 0`. *Proof:* `Psi(c) := c - Hess F(0)^-1 grad F(c)` maps `I^n` to `I^n`, and `c - c' in I^r` (`c, c' in I^n`) implies
`Psi(c) - Psi(c') = -Hess F(0)^-1 [grad F(c) - grad F(c') - Hess F(0)(c - c')] in I^(r+1)`, the bracket being quadratic in `(c, c')` times
`(c - c')`. So the iterates converge and the fixed point is unique. QED (`Hess F(c)` is then invertible, being `Hess F(0)` mod `I`.)
**Lemma T (translation invariance).** Let `F, c, a` be as in s.1.4 and let `p in A^n` be small in one of two senses:
(a) `A = A_0[[q]]`, `F` and `a` defined over `A_0` (or `a` over `A_0[[q]]`), `p in (qA)^n`;  (b) `p in (zA[[z]])^n`.
Define the translated amplitude `b_p(t) := exp( (F(c + t + p) - F(c + t))/z ) a(c + t + p)`, so that
`e^(F(c+t)/z) b_p(t) = e^(F(c+t+p)/z) a(c+t+p)`. Then `b_p` is an admissible amplitude (in `R_n(A_0)[[q]]`, resp. `R_n(A)`), and
`SP^norm_c[F; b_p] = SP^norm_c[F; a]`   (the formal analogue of `int g(x + p) dx = int g(x) dx`).
*Admissibility.* (a) `F(c+t+p) - F(c+t)` has coefficients in `qA`, so its exponential divided by `z` converges q-adically, each
q-coefficient lying in `A_0[[t]][z^-1]`. (b) For `p = z v`: `(F(c+t+zv) - F(c+t))/z = grad F(c+t) . v + O(z)`, and `grad F(c+t) = O(t)`
because `c` is critical; so the exponent lies in the ideal `(t, z)` of `A[[t, z]]` and the exponential converges; `a(c+t+zv)` is in
`A[[t]]((z))`.
*Proof.* Let `s` be a new indeterminate and `B(s) := b_(sp)`. At each fixed `(t-degree, z-power)` (and q-power in case (a)) the
coefficient of `B(s)` is a polynomial in `s` (each `p_i` carries a factor `q`, resp. `z`; in case (b) the factors
`(grad F(c+t) . v)^r` carry `t^r`). By Lemma U every coefficient of `SP^norm_c[F; B(s)]` is a polynomial in `s`, and `d/ds` commutes with
`SP^norm_c[F; .]`. Differentiating `e^(F(c+t)/z) B(s) = e^(F(c+t+sp)/z) a(c+t+sp)` in `s` gives
`d B/ds = Sum_i p_i ( d_i B + z^-1 (d_i F)(c+t) B )`   (`d_i = d/dt_i`; both sides equal `Sum_i p_i` times the `x_i`-derivative of
`e^(F(x+sp)/z) a(x+sp)` divided by `e^(F(c+t)/z)`). By FL1 with `h_i = p_i B(s)` (the `p_i` are constants in `t`; in case (a) FL1 is
applied coefficientwise in `q`, each q-coefficient of `B(s)` lying in `R_n(A_0)[s]`),
`d/ds SP^norm_c[F; B(s)] = 0`. A polynomial over a Q-algebra with zero derivative is constant: `SP^norm[F; B(1)] = SP^norm[F; B(0)]`. QED
**Corollary T' (moving the base point; PROVED).** Let `F_0` over `A_0` have the nondegenerate critical point `c_0`, let
`F = F_0 + delta` with `delta in q A_0[[t]][[q]]`, and `a in A_0[[t]]((z))[[q]]`. Then `F` has a unique critical point `c in c_0 + (qA)^n`
(IFT over `A = A_0[[q]]`), nondegenerate, and in `A_0((z))[[q]]`
`SP^norm_(c_0)[F_0; e^(delta/z) a] = e^((F(c) - F_0(c_0))/z) (Delta(c)/Delta_0)^(-1/2) SP^norm_c[F; a]`,
with `Delta(c) = det(-Hess F(c))`, `Delta_0 = det(-Hess F_0(c_0))`, the binomial series for `(1 + qA)^(-1/2)`, and the exponential series of
`(F(c) - F_0(c_0))/z in q A_0[z^-1][[q]]`. The left side is computed termwise in `q` (s.1.2).
*Proof.* Put `p = c - c_0`. By Lemma T(a) applied to the phase `F_0` and the amplitude `e^(delta/z) a`, the left side equals
`SP^norm_(c_0)[F_0; e^((F(c_0+t+p) - F_0(c_0+t))/z) a(c+t)]`. Expand `F(c + t) = F(c) + t^T H_c t/2 + F_(c,3)(t)` (no linear term) and
`F_0(c_0 + t) = F_0(c_0) + t^T H_0 t/2 + F_(0,3)(t)`. Then `e^(F_(0,3)/z) e^((F(c+t) - F_0(c_0+t))/z) = e^((F(c) - F_0(c_0))/z) e^(t^T (H_c - H_0) t/(2z)) e^(F_(c,3)/z)`:
the `F_(0,3)` cancels. `H_c - H_0` has entries in `qA`, so (G5') turns `E_(-zH_0^-1)[e^(t^T(H_c - H_0)t/(2z)) .]` into
`(det H_c/det H_0)^(-1/2) E_(-zH_c^-1)[.]`, and `E_(-zH_c^-1)[e^(F_(c,3)/z) a(c+t)] = SP^norm_c[F; a]`. Finally
`det H_c/det H_0 = Delta(c)/Delta_0`. QED

## 4. Parameters: Lemma L1c and FL4 (PROVED)

* **(i) Analytic and formal dependence.** Let `F(x; lam)` and `a(x; lam; z)` (coefficientwise in `z`) be holomorphic in `(x, lam)` near
  `(c_0, lam_0)`, `lam in C^p`, with `grad_x F(c_0; lam_0) = 0` and `Hess_x F(c_0; lam_0)` invertible. Then there is a unique holomorphic
  `c(lam)`, `c(lam_0) = c_0`, of critical points (holomorphic IFT), nondegenerate near `lam_0`, and every coefficient of
  `SP^norm_(c(lam))[F(.; lam); a(.; lam)]`, as well as `F(c(lam); lam)` and `Delta(lam)`, is holomorphic near `lam_0`. *Proof:* Lemma U; the
  Taylor coefficients at `c(lam)` are holomorphic in `lam` (Cauchy integrals), `C(lam) = -Hess^-1` is holomorphic. The formal version
  (`A = A_0[[lam - lam_0]]`, critical point from the IFT of s.3.2) is the same statement with "holomorphic" replaced by "in
  `A_0[[lam - lam_0]]`".
* **(ii) Taylor expansion and evaluation.** The Taylor map `O_(lam_0) -> A_0[[lam - lam_0]]` and the evaluation `lam -> lam_0` are ring
  homomorphisms; by (U1) they commute with `SP^norm`. So the Taylor expansion of the analytic SP coefficients is the formal SP over
  `A_0[[lam - lam_0]]` of the Taylor-expanded data, and the value at `lam_0` of each coefficient function is the SP of the data at
  `lam_0`. (Used in chart step (d): "the chart SP at `q = 0` is the value at `q = 0` of its coefficient functions".)
* **(iii) Perturbation in a parameter (T-q).** This is Corollary T': the Taylor expansion in `q` of the SP at the moving critical point
  `c(q)` (with the factors `e^((F(c) - F_0(c_0))/z)` and `(Delta(c)/Delta_0)^(-1/2)`) equals the perturbative expansion at the fixed critical
  point `c_0`, `Sum_d q^d [q^d] SP^norm_(c_0)[F_0; e^(delta/z) a]`. For `delta = q g`, `a` free of `q`:
  `[q^d] = (1/d!) SP^norm_(c_0)[F_0; (g/z)^d a]`. (Used in chart step (e), `g = y^5/prod_(j>=1) x_j`.)
* **(iv) A z-dependent part of the phase (T-z).** Let `g in A[[t]]` with `g(c) = 0` (a constant `g(c)` only multiplies both sides
  by the same symbol `e^(g(c))`). The phase `F + z g` has, after replacing the coefficient `z` by an independent `eps`, the critical
  point `c(eps) in c + (eps A[[eps]])^n`; define the non-perturbative series
  `S(eps; z) := e^((F_eps(c(eps)) - F(c))/z) (Delta_eps/Delta)^(-1/2) SP^norm_(c(eps))[F + eps g; a]`. Then:
  (1) `S(eps; z) = SP^norm_c[F; e^(eps g/z) a]` (Corollary T' with `q = eps`);
  (2) giving `eps` weight 1, every term of `e^(eps g/z)` has weight `>= 1/2` (`g = O(t)`), so by (1) `S` lies in the weight-completed
  ring of series in `eps` and `z` (s.1.2: finitely many terms `eps^d z^k` of each weight `d + k`, weights bounded below), and the
  substitution `eps -> z` is a well-defined, weight-preserving A-linear map;
  (3) `[S(eps; z)]_(eps = z) = SP^norm_c[F; e^g a]`, because the substitution commutes with `E_zC` (it only relabels the coefficients,
  finitely many per weight) and maps `e^(eps g/z)` to `e^g`.
  So "the formal SP with the z-dependent phase `F + zg`" is well defined and equals the SP with `zg` moved into the amplitude as `e^g`.
  This justifies gamma's translation (S4) (GAMMA_PROOFS s.1), the shift `sigma -> sigma + z` of Theorem P(a), and the a-shift in L3.

### 4.5 FL4: the formal SP commutes with d/dq (PROVED)
**Lemma FL4.** Let `F(x; q)`, `a(x; q; z)` be formal in `q - q_0` (or holomorphic in `q`), with the nondegenerate critical point
`c(q)`. Put `S(q) := e^((F(c(q); q) - F(c(q_0); q_0))/z) (Delta(q)/Delta(q_0))^(-1/2) SP^norm_(c(q))[F(.; q); a(.; q)]` (the full SP divided by
its `q`-independent symbols). Then
`dS/dq = e^((F(c(q); q) - F(c(q_0); q_0))/z) (Delta(q)/Delta(q_0))^(-1/2) SP^norm_(c(q))[F; z^-1 (d_q F) a + d_q a]`.
In full-SP language: `d/dq SP[e^(F/z) a dx] = SP[d_q(e^(F/z) a) dx]`; hence also `D_q = q d/dq`, and `P(zD_q) SP[...] = SP[P(sigma) ...]`
when `D_q F = sigma` (the Mellin-Barnes case `F = sigma log(-q) + ...`), for every polynomial `P`.
*Proof.* Write `q = q_0 + eps`. By (iii) at the base point `q_0`, `S(q_0 + eps) = SP^norm_(c(q_0))[F(.; q_0); e^((F(.;q_0+eps) - F(.;q_0))/z) a(.; q_0+eps)]`.
The functional on the right does not depend on `eps`; `d/deps` acts on the amplitude:
`d/deps [e^((F(q_0+eps) - F(q_0))/z) a(q_0+eps)] = e^((..)/z) [z^-1 (d_q F)(q_0+eps) a + (d_q a)(q_0+eps)]`. Apply (iii) again, now to the
amplitude `z^-1 (d_q F) a + d_q a` (an element of `A[[t]]((z))[[eps]]`). Both sides are power series in `eps`; this is the claim at
`q = q_0 + eps`. The holomorphic version follows by Taylor expansion at every `q_0` ((ii)). QED

## 5. Lemma L1a: changes of variables (PROVED)

**Lemma L1a.** Let `F in A[[x - c]]` have the nondegenerate critical point `c`, and let `phi` be a formal map with `phi(c') = c`, i.e.
`phi in (c + (y - c') A[[y - c']])^n`, with `J := D phi(c') in GL_n(A)`. Then `c'` is a nondegenerate critical point of `F o phi`, with
`Hess(F o phi)(c') = J^T H J`, and for every amplitude `a`
`SP^norm_c[F; a] = SP^norm_(c')[F o phi; (a o phi) det D phi / det J]`.
Full SP: `SP_c[e^(F/z) a dx] = SP_(c')[e^((F o phi)/z) (a o phi) det(D phi) dy]` with the branch rule
`det(-J^T H J)^(1/2) := det(J) det(-H)^(1/2)` (for real data with `det J > 0` and positive roots this is the usual change of variables).
Since `A` is arbitrary, this covers analytic maps (via their Taylor series), formal maps and parameter-dependent maps
(`A = A_0[[u]]`, or holomorphic in `u` by s.4 (i)-(ii)), including critical points that move with the parameter.
*Proof.* Translate so that `c = c' = 0`.
*Step 1 (linear maps).* For `phi(y) = J y`: `C' = -(J^T H J)^-1 = J^-1 C J^-T`, `(F o J)_3 = F_3 o J`, and by (G2)
`E_(zC')[e^(F_3(Jy)/z) a(Jy)] = E_(zC)[e^(F_3/z) a]`, the normalised Jacobian being `det J/det J = 1`.
*Step 2 (reduction).* `phi = J psi` with `psi := J^-1 phi`, `D psi(0) = I`, and `det D phi = det J det D psi`. By Step 1 applied to `F`, it
suffices to prove `SP^norm_0[G; b] = SP^norm_0[G o psi; (b o psi) det D psi]` for `G = F o J`, `b = a o J`.
*Step 3 (maps tangent to the identity).* `psi(y) = y + chi(y)`, `chi in (y)^2`. Let `s` be an indeterminate, `psi_s := y + s chi(y)`,
`M_s := D psi_s = I + s D chi` (invertible over `A[s][[y]]`), `m_s := (b o psi_s) det M_s` and `Phi(s) := SP^norm_0[G o psi_s; m_s]`. As
`Hess(G o psi_s)(0) = Hess G(0)` for all `s`, `Phi(s) = E_zC[e^(P_s/z) m_s]` with `P_s := G o psi_s - y^T H y/2` and one fixed `C`. Each
coefficient of `Phi` is a polynomial in `s`, and
`Phi'(s) = E_zC[ e^(P_s/z) ( z^-1 (d_s (G o psi_s)) m_s + d_s m_s ) ]`.
Put `W := M_s^-1 chi` (a vector over `A[s][[y]]`). Two identities in `A[s][[y]]`:
(a) `d_s (G o psi_s) = Sum_k d_k(G o psi_s) W_k`: both sides equal `Sum_i (d_i G)(psi_s) chi_i`, since `Sum_k (M_s)_(ik) W_k = chi_i`;
(b) `d_s m_s = Sum_k d_k(m_s W_k)`: with `adj M = det(M) M^-1`, `m_s W = (b o psi_s) adj(M_s) chi` and
`Sum_k d_k(m_s W_k) = Sum_i (d_i b)(psi_s) (M_s adj(M_s) chi)_i + (b o psi_s) [ tr(adj(M_s) D chi) + Sum_i chi_i Sum_k d_k adj(M_s)_(ki) ]`.
The first term is `det M_s Sum_i (d_i b)(psi_s) chi_i`; `tr(adj(M_s) D chi) = d_s det M_s` (Jacobi's formula, `d_s M_s = D chi`); and the last
sum vanishes by the **Piola identity** `Sum_k d_k adj(D psi)_(ki) = 0`. This is the left side `d_s m_s`.
With (a), (b): `Phi'(s) = SP^norm_0[G o psi_s; Sum_k (d_k h_k + z^-1 d_k(G o psi_s) h_k)]`, `h_k = m_s W_k`, which is 0 by FL1 (phase
`G o psi_s` over `A[s]`). So `Phi(1) = Phi(0)`. QED
*Piola identity (PROVED).* In the formal de Rham complex of `B[[y]]` (any Q-algebra `B`; `d^2 = 0` and the Leibniz rule are formal), the
`(n-1)`-form `omega_i := dpsi_1 ^ .. ^dpsi_i ^ .. ^ dpsi_n` (factor `i` omitted) is closed, because each `dpsi_j` is closed. Its
coefficient on `dy_1..^dy_k..dy_n` is the minor `det M^(i,k)` (row `i`, column `k` deleted), and
`d(f dy_1..^dy_k..dy_n) = (-1)^(k-1) d_k f dy_1..dy_n`. So `Sum_k (-1)^k d_k det M^(i,k) = 0`; multiply by `(-1)^i`:
`Sum_k d_k adj(M)_(ki) = 0`. QED
*Remark (one functional).* The Stage-4 codes compute the same functional: Wick contraction (critic `cr4_chart.py`; gamma's
definition, GAMMA_PROOFS s.1), and the Morse / Lagrange-Buermann form (`rf_sp.py`, `cr4_sp.py`), which is the Wick SP in the Morse
coordinate `y` (phase `F(c) + (H/2) y^2`, a pure Gaussian, so only the amplitude has to be expanded), hence by L1a the Wick SP in the
original coordinate.
### 5.4 FL3 for the SP: the CFU map, termwise in u (PROVED)
In Theorem P(b) (POLY s.6 step 3) and GAMMA_PROOFS s.4-s.5 the map is `Y = X psi(X^5; u)` with `psi in 1 + u Q[X^5][[u]]`, and
`What = f o Y`, `Bhat dX = beta(Y) dY` exactly. At the colour-0 saddle `X = 1` put `t = X - 1`, `Y_c(u) := Y(1; u)` and
`phi(t) := Y(1 + t; u) - Y_c(u) in t Q[[u]][[t]]`; `phi'(0) = Y'(1; u) in 1 + uQ[[u]]` is a unit and the Hessian of the hatted phase
`-f` at `Y_c(u)` is a unit of `Q[[u]]` (at `u = 0`: `f''(1) = 2V(0) + 3 = 5`). So L1a over `A = Q[[u]]` gives
`SP_(X=1)[e^(-What/w) Bhat dX] = SP_(Y=Y_c(u))[e^(-f/w) beta dY]` in `Q[[u]]((w))` (prefactors by the branch rule), i.e. termwise in `u`;
the other colours `X = zeta^-i` likewise over `Q(zeta)[[u]]`. The critical value is the same on both sides, `What(1) = f(Y_c)`.
