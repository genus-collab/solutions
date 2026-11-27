# Stage 5 / chart: the chart argument for (Gamma), written out in full (2026-11-14)

Folder `allgenus/stage5/chart/` (task key `chart`). New files only. AI-written: no human and no referee has read this note.
Labels: **PROVED** (argument written out here, or a quoted theorem whose statement I read, with its number), **VERIFIED** (exact
computation; range and field or primes stated), **CONJECTURAL**, **SKETCH**, **REFUTED**. Finite checks are evidence, never proof.
Notation of Stage 4 (SYN4 = `../../stage4/STAGE4_SYNTHESIS.md`, SYN4b = its part 2): colour 0, `zeta = exp(2 pi i/5)`,
`u = 1 - 5^5 q`, `L = u^(-1/5)`, `b_k = B_2k/(2k(2k-1))`, `h = h_(w_0)`, `w_0 = {1 - zeta^j : j = 1..4} U {-5}`, `h_n = [z^n] h`.
UR, GI, EVID, GAMMA, POLY, QRR as in SYN4. Parts: this file (s.0 statements, s.1 definitions); `CHART_PROOF_2.md` (s.2 the formal
lemmas, with proofs); `CHART_PROOF_3.md` (s.3 the steps (a)-(e), the Z/5 covariance, the proof of Theorem D);
`CHART_PROOF_4.md` (s.4 branches, s.5 what is not proved, s.6 checks, s.7 files). Code and logs: s.6.

## 0. Statements

### 0.1 Theorem Gamma

Put `Lam := C \ {0, 1, zeta, zeta^2, zeta^3, zeta^4}` and, for `L` in `Lam`, `q(L) := (1 - L^-5)/5^5` (so `u = L^-5`, `q != 0`).
`c(z; L) = Sum_n c_n(L) z^n` is the normalised full-Stirling saddle series of the Mellin-Barnes summand at `sigma = L` (s.1.6).

**Theorem Gamma (PROVED in s.3, from the lemmas of s.2; AI-written, unrefereed).**
* (i) Every `c_n` is a rational function of `L` over Q: `c_n in Q[L, L^-1, (1 - L^5)^-1]`.
* (ii) *Chart identity.* For every `L` in `Lam`, `c(z; L) = Nch(z; L)`, the normalised formal stationary-phase series of the chart
  integrand `Omega_ch` (s.3.4) at its critical point `p_ch(L) = (y, x_1..x_4) = (5L, zeta^j - L)`. The coefficients of `Nch` are
  rational functions of `L` over `Q(zeta)` with poles at most at `L = 0, zeta, zeta^2, zeta^3, zeta^4`; they are regular at `L = 1`.
* (iii) **(Gamma).** `c_n` has no pole at `L = 1`, and `c(z; 1) = h(z)`: `c_n(1) = h_n` for every `n >= 0`.
* (iv) *Z/5.* `c_n(zeta L) = zeta^-n c_n(L)`. Hence `c_n` has no pole at any fifth root of unity, `c_n in Q[L, L^-1]`, and
  `c(z; zeta^i) = h(zeta^-i z)`.

### 0.2 Colours: what the Z/5 covariance gives (PROVED, s.3.9)

Colour `i` is the critical point `sigma = zeta^i L`. Let `c^(i)(z; L)` be the normalised series of the same summand at
`sigma = zeta^i L` (the logarithm of `-q` chosen so that this point is critical, s.4). Then
* `c^(i)(z; L) = c(z; zeta^i L) = c(zeta^-i z; L)`, i.e. `c^(i)_n(L) = zeta^(-in) c_n(L)`;
* (Gamma) at colour `i`, `c^(i)(z; 1) = h_(w_i)(z)` with `w_i = {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}`, is EQUIVALENT to (Gamma)
  at colour 0, because `w_i = zeta^i w_0` gives `h_(w_i)(z) = h_(w_0)(zeta^-i z)`.

So the covariance yields nothing beyond colour 0. On the gauge side it is the relation `D_i(z) = D_0(zeta^-i z)` of R0.

### 0.3 Theorem D

**Theorem D (PROVED given the hypotheses below; proof s.3.11).** `D_0 = 1` in `Q[[z]]` to all orders. Hence `D = 1`, and (F_g)
(`D_0 = 1 mod z^(3g-2)`) holds for every `g >= 2`.

Hypotheses (this is the complete list; "Stage 4" means AI-written and AI-refereed in Stage 4, not human-checked):
* **(U)** Theorem U (UR s.3; PROVED there, AI-audited only) and its consequence **R0** (UR s.4): `Q D = P_can Q^m(v(u))`, `D`
  constant diagonal, `D_i(z) = D_0(zeta^-i z)`.
* **(R1)** of UR s.5 (`P_11(0, z) = P_22(0, z)`), used by Lemma N.
* **R2 + R3** (EVID s.2.3-2.4; Stage 4): Lemma N and the connection form, `D = 1 <=> ph^can_n(1) = h_n` for all `n`. Inputs: (U), (R1).
* **R5** Theorem P (POLY s.6 with the poly referee's REFEREE s.3; Stage 4): (a) `Phi^SP` solves HPF and `[w^n] Phi^SP = ph^SP_n(u)`
  is a polynomial; (b) `Phi^SP = Phi^can`. Inputs: (U); GI Theorem G with GI Lemma 1, Lemma V, Lemma Q and (R1); the model rows
  of UR s.1-s.2; the CFU map (POLY s.5, GAMMA_PROOFS s.4); the Brieskorn reduction (POLY s.6 step 2); (FL).
* **(FL)**, the list of the poly referee (C6) plus Stirling's functional equation: formal SP and formal residue kill exact forms;
  both are invariant under formal changes of variables, termwise in `u`; formal SP commutes with `D_q`; `Gamma_St(y+1) = y Gamma_St(y)`
  in `Q[[1/y]]`. The three SP items are PROVED here: (S1), (L1a), (L1c)(iii') of s.2. The residue items and the functional
  equation are not written out here (s.5, item 3).
* **(QRR)** the `q = 0` normalisation of the actual R-matrix, `R(q = 0) = diag(h_(w_i))`; in the flat colour-0 form the QRR WKB
  series has the value `h(z)` at `L = 1` (FORMAL_QUINTIC_genus2 s.0; sources as in SYN4 s.2). It enters through R3.
* **Theorem Gamma** (this note), and the elementary hatted/flat identity `c_n(L) = L^(4n) ph^SP_n(L^-5)` (PROVED in s.3.10).

Theorem D does not use: Lemma P (R1 of the SYN4 table), R4, R6 (Theorem A(i)), any choice of conifold frame, any finite computation.

### 0.4 Consequence

With Theorem D, (F_g) holds for all `g`. By MECHANISM_STAGE4 s.2, Conjecture N(g) then holds for every `g >= 2` modulo (U), the
quoted theorems (Q) of MECHANISM_STAGE4 s.2, (QRR), the unchanged Stage-2 chain, and the Stage-4 inputs listed in s.0.3. No finite
identity remains among the hypotheses.

### 0.5 Status at a glance

| item | label | where |
|---|---|---|
| definitions of SP, `N`, `c(z; L)`; comparison of the Stage-4 definitions | (definitions); comparison PROVED | s.1 |
| (S1), (PG), (L0), (L1a), (L1b), (L1c), (L2) | PROVED (self-contained, algebraic) | s.2 |
| (L3) Stirling = SP of Euler's integral with `x^m` | PROVED modulo DLMF 5.11.1 and DLMF 2.3.8 (read); VERIFIED exactly to `z^40` | s.2.8 |
| steps (a)-(d), Theorem Gamma (i)-(iii) | PROVED | s.3 |
| step (e), coefficientwise agreement with the I-function | PROVED (via (L1c)(iii)); VERIFIED to `q^6 z^10`; not used for Theorems Gamma, D | s.3.7 |
| Z/5 covariance (iv), s.0.2 | PROVED | s.3.9 |
| Theorem D | PROVED given the list of s.0.3 | s.3.11 |
| all of the above | AI-written, unrefereed | - |

## 1. Definitions

### 1.1 Formal integrands

Let `p` be a point of `C^n` and `O_p` the ring of holomorphic germs at `p`. A *formal integrand* at `p` is `e^(Phi/z) A d^n x` with
`Phi` in `O_p` and `A = Sum_(m>=0) A_m z^m` in `O_p[[z]]` (by linearity also `z^-k O_p[[z]]`). The point `p` is a *non-degenerate
critical point* if `dPhi(p) = 0` and `H := Hess Phi(p)` is invertible. With parameters `t` (s.1.4) the data are holomorphic in
`(x, t)`. Everything below is formal in `z`; no integral is ever evaluated, except in the proof of (L3).

### 1.2 The Wick (Feynman) expansion

Write `x = p + xi` and `Phi(p + xi) = Phi(p) + (1/2) xi^T H xi + Phi_3(xi)`, `Phi_3 = O(|xi|^3)`. The *covariance* is `C := -H^-1`.
For a polynomial `P(eta)` the Gaussian expectation `<P>_C` is defined by Isserlis-Wick: `<eta_(i_1) ... eta_(i_2k)>_C` is the sum
over the perfect matchings of `{1..2k}` of the products of the `C_(i_a i_b)`; odd moments are 0. Equivalently
`<P>_C = [exp((1/2) Sum_(i,j) C_ij d_i d_j) P](0)`, or the recursion `<eta_i P>_C = Sum_j C_ij <d_j P>_C`, `<1> = 1`.
* **Coefficients.** `N_n[Phi, A; p] := Sum_(m=0..n) Sum_(V=0..2(n-m)) (1/V!) < [Phi_3^V A_m]_(2(n-m+V)) >_C`, where `[F]_d` is the
  homogeneous part of degree `d` of `F(xi)`, read as a polynomial in `eta`. Equivalently: put `xi = z^(1/2) eta` in
  `exp(Phi_3(xi)/z) A(p + xi, z)`, expand in powers of `z^(1/2)` (the coefficient of `z^(r/2)` is a polynomial in `eta` of parity
  `r`), apply `< >_C` termwise; the half-integer powers drop out. The sum is finite: `[Phi_3^V]_d = 0` for `d < 3V`.
* `N[Phi, A; p](z) := Sum_n N_n z^n`; `N_0 = A_0(p)`. The **normalised series** is `Nhat := N/N_0` (when `A_0(p) != 0`).
* **The SP symbol.** `SP_p[e^(Phi/z) A d^n x] := (2 pi z)^(n/2) det(-H)^(-1/2) e^(Phi(p)/z) N[Phi, A; p](z)`: Gaussian factor
  `(2 pi z)^(n/2) det(-H)^(-1/2)` (a branch of the square root is part of the datum), critical value `Phi(p)`, series `N`.

Facts used below (all immediate from the definition):
* **(F1)** `N_n` is a polynomial with rational coefficients in the Taylor coefficients of `Phi` at `p` of orders `3..2n+2`, those of
  `A_m` of orders `<= 2(n-m)`, and the entries of `C`. (A term with `V` vertices of degrees `k_i` and an amplitude monomial of
  degree `j` has `Sum k_i + j = 2(n - m + V)`, so `k_i <= 2n + 2` and `j <= 2(n - m)`.)
* **(F2)** `N` does not see the constant `Phi(p)`: adding a constant to `Phi` changes only the critical value.
* **(F3) Homogeneity.** For `mu` in `C^x`: `N[mu Phi, A(., z)](z) = N[Phi, A(., mu z')](z')` with `z' = z/mu`. (The integrand
  `e^(mu Phi/z) A` depends on `(Phi, z)` through `Phi/z'`; with `eta' = mu^(1/2) eta` the covariance `-(mu H)^-1` becomes `-H^-1`.)
* **(F4) Analytic meaning (remark only; used in no proof).** For real data with `H` negative definite, `SP_p` is the asymptotic
  expansion of `Int e^(Phi/z) A dx` over a small ball as `z -> 0+` (Laplace's method; e.g. Hoermander, ALPDO I, Thm 7.7.5, whose
  operators `L_j` are exactly this Wick formula; statement recalled from memory, not re-read).

### 1.3 Rings

The same formulas make sense with coefficients in any commutative Q-algebra `R` in which `H` is invertible (`R = O_(t_0)` for
parameters; `R = Q[[w]]/(w^D)` or `Q(zeta)[[q]]` in the checks). All lemmas of s.2 hold over such `R` with the same proofs.

### 1.4 Parameters

Data depending holomorphically on `t` near `t_0`, with a non-degenerate critical point `p_0` at `t_0`, have a unique holomorphic
family of critical points `p(t)`, `p(t_0) = p_0` (implicit function theorem, since `det H != 0`), non-degenerate for `t` near `t_0`.

### 1.5 Exact forms

`SP_p[d(e^(Phi/z) g)]` means `SP_p[Sum_i d_i(e^(Phi/z) g_i) d^n x]` for an `(n-1)`-form `g = Sum_i (-1)^(i-1) g_i dx_1..^i..dx_n`
with `g_i` in `O_p[[z]]`; its amplitude is `z^-1 Sum_i ((d_i Phi) g_i + z d_i g_i)`.

### 1.6 The series c(z; L), exactly

Fix `L` in `Lam` and an open disc `Delta` around `L` containing neither 0 nor any `zeta^j`. On `Delta` fix holomorphic branches of
`log(5 sigma)` and `log(zeta^j - sigma)`, `j = 0..4`, and put
`lambda := -5 log(5L) + Sum_j log(zeta^j - L)`, so that `e^lambda = prod_j (zeta^j - L)/(5L)^5 = (1 - L^5)/(5^5 L^5) = -q(L)`:
`lambda` is the logarithm of `-q` that makes `sigma = L` critical (s.4, (B1)). Define on `Delta`
* `S_0(sigma) := lambda sigma + 5 sigma log(5 sigma) - 5 sigma + Sum_j [(zeta^j - sigma) log(zeta^j - sigma) - (zeta^j - sigma)]`,
* `S_1(sigma) := (1/2) log(5 sigma) - (1/2) Sum_j log(zeta^j - sigma)`,
* `S_2k(sigma) := b_k [(5 sigma)^(1-2k) + Sum_j (zeta^j - sigma)^(1-2k)]`, `k >= 1`.

Then `S_0' = lambda + 5 log(5 sigma) - Sum_j log(zeta^j - sigma)`, `S_0'(L) = 0`, and
`S_0'' = 5/sigma + Sum_j 1/(zeta^j - sigma) = 5/(sigma (1 - sigma^5))`, which is non-zero at `L`.

**Definition.** `c(z; L) := Nhat[S_0, A_St; L](z)`, one variable `sigma`, measure `d sigma`, amplitude
`A_St(sigma, z) := exp(S_1(sigma)) exp(Sum_(k>=1) z^(2k-1) S_2k(sigma))` (in `O_L[[z]]`).

*Why this is the summand with all six Gamma factors Stirling-expanded.* For `a` in `{5 sigma, zeta^j - sigma}` read `log(a/z)` as
`log a - log z`, and `Gamma_St(1 + w) := w Gamma_St(w)` (equivalently, expand Stirling's formula at `1 + w` in powers of `1/w`;
the two agree by Stirling's functional equation). Then, for each `sigma` in `Delta`,

    e^(lambda sigma/z) Gamma_St(1 + 5 sigma/z) prod_(j=0..4) Gamma_St((zeta^j - sigma)/z) = (2 pi)^3 z^2 exp(S_0/z + S_1 + Sum_k z^(2k-1) S_2k),

because the coefficient of `(log z)/z` is `-[5 sigma + Sum_j (zeta^j - sigma)] = -Sum_j zeta^j = 0` and that of `log z` is
`-1/2 + 5/2 = 2` (both identities in `sigma`; VERIFIED symbolically by `ch_zpow.py` (i)-(iii)). So `c(z; L)` is the normalised formal
SP at `sigma = L` of `Gamma(1+5s) prod_j Gamma(zeta^j/z - s) (-q)^s d sigma` with all six Gamma factors replaced by `Gamma_St` and
`(-q)^s := e^(lambda sigma/z)`.

*Taylor data (used in s.3 and in the code).* With `P_k(L) := Sum_j (zeta^j - L)^-k` (symmetric in the `zeta^j`, so in `Q(L)`;
generating function `Sum_k P_k t^(k-1) = 5 (L+t)^4/(1 - (L+t)^5)`):
* `S_0(L+t) - S_0(L) = Sum_(k>=1) t^(k+1) [5 (-1)^(k+1) L^-k + P_k]/(k(k+1))`;
* `S_1(L+t) - S_1(L) = Sum_(k>=1) t^k [(-1)^(k+1) L^-k + P_k]/(2k)`;
* `[t^r] S_2k(L+t) = b_k binom(1-2k, r) [5^(1-2k) L^(1-2k-r) + (-1)^r P_(2k-1+r)]`.

By (F1), `c_n(L)` is a polynomial over Q in these numbers and in `-1/S_0''(L) = -L(1 - L^5)/5`. Since
`P_k in (1 - L^5)^-k Q[L]`, this proves **Theorem Gamma (i)**: `c_n in Q[L, L^-1, (1 - L^5)^-1]`. The value does not depend on
the choices (disc, branches): the Taylor coefficients of order `>= 1` of `S_0`, `S_1`, `S_2k` do not (s.4).

### 1.7 The Stage-4 definitions are this object (q versus -q, dx/x versus dx)

* **GAMMA s.2** (and `gm_fs.py`): the same integrand `Gamma(5s+1) prod_j Gamma(zeta^j/z - s) (-q)^s`, the same `S_0'`, `S_1`,
  `S_2k`, measure `d sigma`, the same normalisation `c_0 = 1`. Identical to s.1.6.
* **POLY s.5** (and the critic's `cr4_sp.py`, which implements it at fixed `u0`): the I-function summand
  `q^s Gamma(1+5s)/prod_j Gamma(1 + s - zeta^j/z)`, then `x = 1/sigma`, `x = sX` (`s = 1/L`) and `z = s^4 w`. Three differences,
  none of which changes the normalised series:
  * *Reflection.* For `y = (zeta^j - sigma)/z`, Stirling's formula gives, as formal expressions,
    `log Gamma_St(y) + log Gamma_St(1 - y) = log 2 pi + (y - 1/2)(log y - log(-y))` (the `b_k` terms cancel: `1 - 2k` is odd), and
    `log y - log(-y) = i pi eps_j`, `eps_j` odd. So `prod_j 1/Gamma_St(1 + (sigma - zeta^j)/z)` is a constant times
    `e^(c/z) e^(i pi (Sum_j eps_j) sigma/z) prod_j Gamma_St((zeta^j - sigma)/z)`. Together with `q^s` versus `(-q)^s`
    (`log q = log(-q) + i pi eps`), POLY's integrand is a constant times GAMMA's with `lambda` replaced by `lambda + i pi e`,
    `e = eps + Sum_j eps_j` even. Requiring `L` to be critical forces `e = 0` (the branches of POLY's logarithms are chosen that
    way). The normalised series depends only on the Taylor coefficients of order `>= 2` of the phase and on the normalised amplitude.
  * *Variables.* `sigma = L/X` is an analytic change of variables near `X = 1`: (L1a).
  * *Rescaling* `z = L^-4 w`: (F3). It gives `c_n(L) = L^(4n) ph^SP_n(L^-5)` (s.3.10). POLY divides by `sqrt(2 pi w/5) e^(-What(1)/w)`,
    which is the Gaussian factor times the amplitude value `(1-u)^(-1/2)`: `Phi^SP` is normalised.
* **`independent_check/gamma_interp.py`**: the sigma-form in POLY's sign convention, amplitude
  `sqrt(5 sigma/(sigma^5 - 1)) exp(Sum_k b_k z^(2k-1) [(5 sigma)^(1-2k) - Sum_j (sigma - zeta^j)^(1-2k)])`, which is `A_St` up to a
  constant (odd exponents), measure `d sigma`: the same `c(z; L)`.
* **`critic/cr4_chart.py`**: the chart integrand of s.3.4 with amplitude `prod_(j>=1) 1/x_j` and measure `dy`; it drops the constant
  `2 pi kappa (-q)^(1/z)` (s.3.4) and compares with the QRR WKB `R_n(L) = L^(4n) ph_n(L^-5)`, not with `c_n`; they agree by (Gamma)
  and Theorem P(a) (both are colour-0 WKB series with polynomial coefficients, so `c/R` is a series in `z` alone, equal to 1 at
  `L = 1`). Its "Jacobian `|x_0|`" is the factor `kappa X_0`, `kappa = +-i`, of (L2). SYN4 s.1.4 writes the chart prefactor as
  `2 pi z (-q)^(1/z)`; with `Omega_7 = z^-1 e^(Phi/z) ...` it is `2 pi kappa (-q)^(1/z)`. A power of `z` and a constant: harmless.
* **q versus -q.** In every form the exponentiated critical equation is `-q (5 sigma)^5 = prod_j (zeta^j - sigma)`, i.e.
  `u sigma^5 = 1`. With `+q` the critical points would satisfy `sigma^5 (1 + 5^5 q) = 1` (control in `ch_hess.py`).
* **dx/x versus dx.** `Gamma(1 + 5s) = Int y^(5s) e^-y dy` carries `dy`; `Gamma(w) = Int x^w e^-x dx/x` carries `dx/x`. Replacing
  `dy` by `dy/y` means `Gamma(5s)` for `Gamma(1 + 5s)`; it changes `c` from `z^2` on, and the consistently changed pair agrees
  (`ch_chart.py` controls). (L3) is VERIFIED for both measures.
* **My code** (`ch_phases.sigma_data`) uses the GAMMA sigma-form. It gives `c_1 = -(3/20) L^4`, `c_2 = (9/800) L^8` and reproduces
  the stored value `c_5(2) = 11495493197/31500` (GAMMA logs, SYN4 s.3 E3): `ch_sigma.py`.
