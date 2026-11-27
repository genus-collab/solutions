# DIGEST_D2: Theorem Gamma, LMI, formal stationary phase, free coefficients, generalisation, (QRR) sources

Design digest for the note `allgenus/note/` (key read_D2). Six pieces: this file (s.0 sources, s.1 definitions);
`DIGEST_D2_2.md` (s.2 statements, part A: Theorem Gamma, LMI and consequences); `DIGEST_D2_3.md` (s.2 part B: formal
stationary-phase lemmas for the appendix; (QRR)); `DIGEST_D2_4.md` (s.3 proofs of part A); `DIGEST_D2_5.md` (s.3 proofs of
part B; s.4 conventions); `DIGEST_D2_6.md` (s.5 evidence, s.6 literature, s.7 open issues). No LaTeX is written here.

**Abbreviations** (paths relative to `HVQ/allgenus/`, HVQ = `formalquintic/hodge_vertex_quintic`):
M5 = `MECHANISM_STAGE5.md`; M4 = `MECHANISM_STAGE4.md`; SYN4 = `stage4/STAGE4_SYNTHESIS.md`;
CP, CP2, CP3, CP4 = `stage5/chart/CHART_PROOF.md`, `_2`, `_3`, `_4` (the 7-variable chart, lemmas (S1)..(L3));
RF, RF2 = `stage5/chart/referee_formal/REFEREE.md`, `_2`; RA, RA2 = `stage5/chart/referee_adversary/REFEREE.md`, `_2`
(RA s.3 = Theorem LMI and (L2')); LF, LF2 = `stage5/referees/lmi_formal/REFEREE.md`, `_2`; LA, LA2 =
`stage5/referees/lmi_adversary/REFEREE.md`, `_2`; FSL, FSL2..FSL5 = `stage5/lemmas/FORMAL_SP_LEMMAS.md`, `_2`.. `_5`;
RL, RL2 = `stage5/lemmas/referee/REFEREE.md`, `_2`; GEN, GEN2 = `stage5/general/GENERAL.md`, `_2`; LIT, LIT2, LIT3 =
`stage5/lit/LIT.md`, `_2`, `_3`; SYN5, SYN5b = `stage5/STAGE5_SYNTHESIS.md`, `_2`; RS, RS2, RS3 =
`stage5/referees/REFEREES_SYNTHESIS.md`, `_2`, `_3`; IC = `stage5/independent_check/README.md`.
Correction tags as in RS s.0: RF-C1..C6 (RF2 s.6), RL-C1..C8 (RL2 s.5), LF-C1..C6 (LF s.0), LA-C1..C7 (LA2 s.6);
adjudications RS A1-A10, SYN5b K1-K14 and O1-O19; RS2 items 1-48 (the merged list of corrections).
**Status vocabulary** (as in the sources): PROVED (argument written out), with the referee count R2 / R1 / R0 (refereed twice /
once / not at all, always by AI agents; no human has checked anything); VERIFIED (finite exact computation); SKETCH; CONJECTURAL;
quoted (literature, with the statement number and how it was read).

## 0. Sources read, and what was not read

Read in full (every section):
* M5 (all, s.1-s.5); M4 (all, s.1-s.5).
* CP s.0-s.1; CP2 s.2.1-s.2.8; CP3 s.3.1-s.3.11; CP4 s.4-s.7.
* RF s.0-s.3; RF2 s.4-s.9. RA s.0-s.3; RA2 s.4-s.7.
* LF s.0-s.4; LF2 s.5-s.7. LA s.0-s.3; LA2 s.4-s.8.
* FSL s.0-s.2; FSL2 s.3-s.5; FSL3 s.6-s.8; FSL4 s.9-s.10; FSL5 s.11-s.12. RL s.0-s.3; RL2 s.4-s.9.
* GEN s.0-s.4; GEN2 s.5-s.9. LIT s.0-s.1; LIT2 s.2-s.3; LIT3 s.4-s.6.
* SYN5 s.0-s.4; SYN5b s.5-s.10. RS s.0-s.4; RS2 s.5-s.6; RS3 s.7-s.9. IC (README).
Read in part (context only): SYN4 s.1.1-s.1.4 (first lines of s.1.4); `stage2/STAGE2_SYNTHESIS_2.md` s.6.3;
`stage4/STAGE4_SYNTHESIS_2.md` s.6.2-s.6.4; the template note `resolvedconifold/hodge_vertex_conifold/note/`: OUTLINE.md s.0
(binding rules, s.0.1-s.0.3) and LABELS.md (head and the s01-s03 tables); the root `refs.bib`, searched (read-only) for keys.
**Not read** (anything there reaches this digest only as quoted by the files above):
* every script and log (I used the reports' descriptions of them; numbers in s.5 are as the reports state them);
* `stage5/pb/PB_AUDIT*.md`, `stage5/critic/CYC_PROOF.md`, `stage5/referees/cyc_*`, `stage5/referees/pb_steps/` (Theorem P(b),
  the glue, (CYC): other digests);
* Stage-4 worker notes (GAMMA, GAMMA_PROOFS, POLY, EVID), UR (`natural_form/theory/UNIFORM_REDUCTION.md`), GI, NATURAL_FORM,
  FORMAL_QUINTIC_genus2, MECHANISM_STAGE2, MECHANISM, SYNTHESIS_STAGE01;
* the cited papers and DLMF themselves (no web access was used by me; every statement number in s.6 is as LIT, CP, RF, RA, LF,
  LA, FSL, RL or SYN5 report having read it, with the mode of reading they report).

## 1. Definitions (final form; names are proposals; notation to be fixed in NOTATION.md)

**D1. Basic parameters.** `zeta = exp(2 pi i/5)`; `u = 1 - 5^5 q`; `L` with `u L^5 = 1` (`L = u^(-1/5)`, colour 0 near `q = 0`
means `L` near 1). The independent variable is `L`; put `q(L) := (1 - L^-5)/5^5` (single-valued, no fifth root taken; CP4 s.4
(B9)) and `Lam := C \ {0, 1, zeta, zeta^2, zeta^3, zeta^4}` (CP s.0.1). `Lam1(L) := prod_(j=1..4) (zeta^j - L) = 1 + L + L^2 + L^3 + L^4`,
so `prod_(j=0..4) (zeta^j - L) = 1 - L^5 = (1 - L) Lam1(L)` (CP3 header). Sources: CP s.0.1, CP3 header, CP4 s.4.
*Naming clash:* `Lam` (this set), `Lam(z)` (GEN s.3.3), `Lam = -log(1 - u)` (FSL4 s.10.4), `Lam(Y)` (FSL4 s.9.3); choose distinct
symbols in NOTATION.md.

**D2. Bernoulli data and Stirling forms.**
* `b_k := B_2k/(2k(2k-1))`, `k >= 1` (CP s.0, LIT s.0); `b_1 = 1/12`.
* `beta_n(theta) := (-1)^(n+1) B_(n+1)(theta)/(n(n+1))`, `n >= 1`, `B_m(theta)` the Bernoulli polynomials; `beta_(2k-1)(0) = b_k`,
  `beta_(2k)(0) = 0` (GEN s.1; LA s.2.3; LF s.4.2). The sources call the shift parameter `h`: rename it (clash with `h = h_(w_0)`).
* The formal Stirling form `Gamma_St(Y) := exp((Y - 1/2) log Y - Y + (1/2) log 2pi + Sum_(k>=1) b_k Y^(1-2k))` (FSL4 s.9.1 (iii));
  with argument `a/z` one reads `log(a/z) := log a - log z` (CP s.1.6; CP4 (B7)), and `Gamma_St(1 + w) := w Gamma_St(w)`
  (CP s.1.6; consistent with Stirling at `1 + w` by the formal functional equation, FSL4 s.9.1 (iv)).
* `Poch_M(e) := prod_(i=0..M-1) (1 + i e)` for `M >= 0`, `:= prod_(i=1..-M) (1 - i e)^-1` for `M < 0` (CP2 s.2.8).
* `St_m(w) := exp((1/w + m - 1/2) log(1 + m w) - m + Sum_k b_k w^(2k-1) (1 + m w)^(1-2k))` in `Q[[w]]` (FSL4 s.9.1);
  `St_m = Poch_m St_0` and `St_0(w) = exp(Sum_k b_k w^(2k-1))` (RL s.3.1 table; FSL4 s.9.3 A4).

**D3. QRR weights and the value series.** For a finite multiset `w` of non-zero numbers,
`h_w(z) := exp(-Sum_(k>=1) b_k z^(2k-1) Sum_(a in w) a^(1-2k))` (LIT s.0; M4 s.1). Colour weights
`w_i := {zeta^i - zeta^j : j != i} U {-5 zeta^i}`; `h := h_(w_0)`, `w_0 = {1 - zeta^j : j = 1..4} U {-5}`; `h_n := [z^n] h`.
Equivalent form: `T_k := 5^(1-2k) + Sum_(j=1..4) (zeta^j - 1)^(1-2k) = -Sum_(a in w_0) a^(1-2k)` (odd exponent), so
`h(z) = exp(Sum_k b_k T_k z^(2k-1))` (CP3 s.3.6). `w_i = zeta^i w_0` gives `h_(w_i)(z) = h(zeta^-i z)` (CP s.0.2, CP3 s.3.9).
Values (RA2 s.4.2 C0; LF2 Q1; re-computed here from the definition): `h_1..h_5 = -3/20, 9/800, -269/80000, 2823/6400000,
-257071/8064000000`. With the weight `+5` instead of `-5`: `[z] log h = -11/60` (CP3 s.3.6).

**D4. The orbifold-point value and the free coefficients.**
`Ginf(z) := exp(Sum_(l>=1) (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`, "the Stirling factor of `prod_j Gamma(zeta^j/z + 1/5)`"
(SYN4 s.1.1); `[z^5] Ginf = B_6(1/5)/6 = 4537/3937500` (SYN4 s.1.1; LF2 Q4; re-computed here). Free coefficient (SYN4 s.1.1):
`f_m := [u^(4m)] ph_(5m)`, the top coefficient of `ph_(5m)` (`deg ph_n <= floor(4n/5)`, SYN4 s.1.1); `f^SP_m := [u^(4m)] ph^SP_(5m)`
(RA s.3.4), `f^can`, `f^QRR` likewise for the canonical and QRR columns (SYN5b K7). `ph^SP_n` see D10.

**D5. Formal integrands and the Wick functional** (CP s.1.1-s.1.2; same object as FSL s.1.3-s.1.4, Def. 1.2 = CP's `N`, RL-C7).
* A formal integrand at `p in C^n` is `e^(Phi/z) A d^n x` with `Phi` holomorphic at `p` and `A = Sum_(m>=0) A_m z^m`, `A_m`
  holomorphic at `p` (by linearity also `z^-k O_p[[z]]`). `p` is a non-degenerate critical point if `dPhi(p) = 0` and
  `H := Hess Phi(p)` is invertible. Covariance `C := -H^-1`.
* Gaussian expectation `<P>_C` (Isserlis-Wick: sum over perfect matchings of products of `C_ij`; equivalently
  `[exp((1/2) Sum C_ij d_i d_j) P](0)`, or the recursion `<eta_i P> = Sum_j C_ij <d_j P>`, `<1> = 1`).
* `N_n[Phi, A; p] := Sum_(m=0..n) Sum_(V=0..2(n-m)) (1/V!) < [Phi_3^V A_m]_(2(n-m+V)) >_C`, `Phi_3 := Phi(p + xi) - Phi(p) - (1/2) xi^T H xi`,
  `[F]_d` the homogeneous part of degree `d`. Equivalently: `xi = z^(1/2) eta`, expand `exp(Phi_3(xi)/z) A(p + xi, z)` in powers of
  `z^(1/2)`, apply `< >_C` termwise; half-integer powers drop out. `N := Sum N_n z^n`, `N_0 = A_0(p)`.
* Normalised series `Nhat := N/N_0` (requires `A_0(p) != 0`). Symbol (bookkeeping only):
  `SP_p[e^(Phi/z) A d^n x] := (2 pi z)^(n/2) det(-H)^(-1/2) e^(Phi(p)/z) N` (the branch of the root is part of the datum).
* Coefficient rings (CP s.1.3): any commutative Q-algebra in which `H` is invertible; FSL s.1.2 gives the weighted completed ring
  `R_n(A)` (`wt(t^alpha z^j) = |alpha|/2 + j`, finitely many terms of weight `<= W`; Lemma 1.1) and the functional `E_zC` (s.1.3).
* Parameters (CP s.1.4): data holomorphic in `(x, t)` with a non-degenerate critical point at `t_0` have a unique holomorphic family
  of critical points `p(t)` (implicit function theorem), non-degenerate near `t_0`.
* Exact forms (CP s.1.5): amplitude `z^-1 Sum_i ((d_i Phi) g_i + z d_i g_i)`.

**D6. The Mellin-Barnes summand and the series `c(z; L)`** (CP s.1.6, the definition of record; = GAMMA s.2 = POLY s.5 = the
Stage-4 objects, CP s.1.7, RF s.1 CONFIRMED). Summand (variable `sigma = z s`):
`Gamma(1 + 5s) prod_(j=0..4) Gamma(zeta^j/z - s) (-q)^s d sigma`, all six Gamma factors replaced by `Gamma_St`, `(-q)^s := e^(lambda sigma/z)`.
Fix `L in Lam`, a disc `Delta` around `L` avoiding `0` and the `zeta^j`, branches of `log(5 sigma)`, `log(zeta^j - sigma)` on `Delta`,
and `lambda := -5 log(5L) + Sum_j log(zeta^j - L)` (so `e^lambda = (1 - L^5)/(5^5 L^5) = -q(L)`; this is the logarithm of `-q` that makes
`sigma = L` critical, CP4 (B1)). Put
* `S_0(sigma) := lambda sigma + 5 sigma log(5 sigma) - 5 sigma + Sum_j [(zeta^j - sigma) log(zeta^j - sigma) - (zeta^j - sigma)]`;
* `S_1(sigma) := (1/2) log(5 sigma) - (1/2) Sum_j log(zeta^j - sigma)`;
* `S_2k(sigma) := b_k [(5 sigma)^(1-2k) + Sum_j (zeta^j - sigma)^(1-2k)]`, `k >= 1`;
* `A_St(sigma, z) := exp(S_1(sigma)) exp(Sum_(k>=1) z^(2k-1) S_2k(sigma))`.
Then `S_0'(L) = 0`, `S_0'' = 5/sigma + Sum_j 1/(zeta^j - sigma) = 5/(sigma (1 - sigma^5))`, and **`c(z; L) := Nhat[S_0, A_St; L](z)`**, one
variable, measure `d sigma`, `c_0 = 1`. Justification (CP s.1.6; RF s.1 CONFIRMED): with `log(a/z) = log a - log z`,
`e^(lambda sigma/z) Gamma_St(1 + 5 sigma/z) prod_j Gamma_St((zeta^j - sigma)/z) = (2 pi)^3 z^2 exp(S_0/z + S_1 + Sum_k z^(2k-1) S_2k)`, since the
coefficient of `(log z)/z` is `-Sum_j zeta^j = 0` and that of `log z` is `-1/2 + 5/2 = 2`, identically in `sigma`.
*Taylor data* (CP s.1.6), with `P_k(L) := Sum_j (zeta^j - L)^-k`, `Sum_k P_k t^(k-1) = 5 (L + t)^4/(1 - (L + t)^5)`, `P_k in (1 - L^5)^-k Q[L]`:
`S_0(L+t) - S_0(L) = Sum_(k>=1) t^(k+1) [5 (-1)^(k+1) L^-k + P_k]/(k(k+1))`; `S_1(L+t) - S_1(L) = Sum_(k>=1) t^k [(-1)^(k+1) L^-k + P_k]/(2k)`;
`[t^r] S_2k(L+t) = b_k binom(1-2k, r) [5^(1-2k) L^(1-2k-r) + (-1)^r P_(2k-1+r)]`; `-1/S_0''(L) = -L(1 - L^5)/5`.
*Colour i* (CP s.0.2): `c^(i)(z; L)` := the normalised series of the same summand at `sigma = zeta^i L` (the logarithm of `-q` chosen so
that this point is critical).
*Other forms of the same series* (CP s.1.7): POLY's I-function summand `q^s Gamma(1 + 5s)/prod_j Gamma(1 + s - zeta^j/z)` (reflection,
or RF-C5's direct argument, s.3 of DIGEST_D2_4); the hatted variables `x = 1/sigma = X/L`, `z = L^-4 w` (D10).
*Motivation (GEN s.2 Lemma 2, PROVED as a residue identity, SYN5b O3; R0):* for the summand
`a(s) = prod_j Gamma(zeta^j/z - s) Gamma(1 + 5s) (-q)^s`, `Res_(s = zeta^i/z + d) a / Res_(s = zeta^i/z) a = q^d I_d(p_i)`, the degree-`d`
coefficient of the I-function at the fixed point `p_i`. (GEN's general form: `Q = (-1)^(N + l_conc) z^(Sum c_f) q`; here `Q = -q`.)

**D7. LMI data** (RA s.3.1; LF s.1, LF-C3; LA s.1.1, LA-C3). Fix `s* in C`, `alpha0 in C^x`, `M0 in Z`, `g` holomorphic near
`s*` with `g_k := g^(k)(s*)/k!`, `F in O_(s*)[[z]]` (every `F_m` holomorphic near `s*`) with `F(s*, 0) != 0`;
`a0(sigma) := alpha0 (sigma - s*)` (EXACTLY linear). Let `U*` be a punctured disc around `s*` on which `g` and the `F_m` are
holomorphic and `alpha0 + (L - s*) g''(L) != 0` (equivalently `S''(L) != 0`). For `L in U*` choose a branch of `log a0` near `L`:
`lambda(L) := -g'(L) - alpha0 log a0(L)`; `S(sigma) := lambda sigma + g(sigma) + a0 log a0 - a0`;
`A(sigma, z) := a0^(M0 - 1/2) Poch_M0(z/a0) exp(Sum_k b_k (z/a0)^(2k-1)) F(sigma, z)` (the Stirling form of `Gamma(M0 + a0/z)` times
`F`, stripped of `(2 pi z)^(1/2)` and of the phase, LF s.1); **`c^F(z; L) := Nhat[S, A; L]`**. Also `delta := -M0 z/alpha0`,
`gt(x) := g(s* + x) - g(s*) - g'(s*) x`, chart parameter `Q := e^(-lambda/alpha0) = a0(L) e^(g'(L)/alpha0)` (single-valued,
LF s.3.2), `mu_* := -g'(s*)/alpha0`. For `M0 in C` (LA s.2.3): `A := a0^(M0 - 1/2) exp(Sum_n beta_n(M0) (z/a0)^n) F`.
*Naming clashes:* `Q` (here) versus GEN's `Q`; `delta` (here) versus `delta` in (PG); `h` (LF/LA shift) versus `h(z)`.

**D8. The 7-variable integrand and the chart** (CP3 s.3.1, s.3.4). `Omega_7 := z^-1 e^(Phi/z) (x_0 x_1 x_2 x_3 x_4)^-1 d sigma dy dx_0..dx_4`,
`Phi := lambda sigma - y + 5 sigma log y + Sum_(j=0..4) (-x_j + (zeta^j - sigma) log x_j)`, critical point `p(L) = (L, 5L, zeta^j - L)`.
`Nhat_7(z; L)` := its normalised series. `G := lambda + 5 log y - Sum_(j=0..4) log x_j`, `f := -y + Sum_j (-x_j + zeta^j log x_j)`
(`Phi = sigma G + f`); `X_0(w) := -q y^5/prod_(j>=1) x_j`, `w = (y, x_1..x_4)`. Chart phase
`psi(y, x_1..x_4) := 5 log y - y + Sum_(j=1..4) (-x_j + (zeta^j - 1) log x_j) + q y^5/(x_1 x_2 x_3 x_4)`; chart integrand
`Omega_ch := e^(psi/z) dy prod_(j=1..4) dx_j/x_j`; critical point `p_ch(L) = (5L, zeta^1 - L, .., zeta^4 - L)`, where
`K := q y^5/prod x_j = L - 1`; `Nch(z; L)` := its normalised series.

**D9. I-function coefficients** (CP3 s.3.7): `I_d(z) := (1/d!) prod_(k=1..5d) (5/z + k) / prod_(j=1..4) prod_(i=1..d) ((1 - zeta^j)/z + i)`,
the `d`-th coefficient of the colour-0 I-function `z q^(1/z) Sum_d q^d prod_(k<=5d) (5 + kz)/prod_(j=0..4) prod_(k<=d) (1 - zeta^j + kz)`
(FORMAL_QUINTIC_genus2 s.1 as quoted by CP3 s.3.7; the `j = 0` factor gives `d! z^d`, RF s.3.7).

**D10. The hatted series** (POLY s.5 as quoted by CP s.1.7, CP3 s.3.10, RF s.1): `Phi^SP(w; u) = Sum_n w^n ph^SP_n(u)` is the normalised
SP series of POLY's hatted integrand (variables `x = 1/sigma = X/L`, `z = L^-4 w`, saddle `X = 1`). RF s.1 re-derived its data: phase
derivative `d/dX [L^4 S(L/X)] = (1/(u X^2)) log((1 - u X^5)/(1 - u))`, amplitude including the Jacobian
`(1 - u X^5)^(-1/2) exp(Sum_k b_k (uXw)^(2k-1) [5^(1-2k) - 5 Sum_P C(2k-2+5P, 5P) u^P X^(5P)])`. The full definition belongs to the
Theorem-P digest; here only the identity of D2_2 GAM-G is used.
