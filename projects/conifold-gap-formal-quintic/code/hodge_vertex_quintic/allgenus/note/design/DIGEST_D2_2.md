# DIGEST_D2, part 2: s.2 Statements, part A (Theorem Gamma, LMI, consequences)

Continuation of `DIGEST_D2.md` (abbreviations, s.0, s.1 definitions D1-D10). Proposed label prefix `gam` for the section
"the value at q = 0" (Theorem Gamma); the appendix lemmas (prefix `fsp`) are in `DIGEST_D2_3.md`. Each entry: label, type;
statement in final corrected form; hypotheses and scope; source; status; corrections applied; dependencies.

## 2A.1 The series and its forms

**GAM-1** `prop:gam-rational` (prop). *Statement.* For `L in Lam`, `c(z; L)` (D6) does not depend on the disc `Delta` nor on the
branches of the logarithms; every `c_n` lies in `Q[L, L^-1, (1 - L^5)^-1]`. More precisely `c(z; L)` depends only on the Taylor
coefficients of order `>= 2` of the phase at `L` and on the normalised amplitude, so `lambda` and the sign convention for `q` never
enter it; they only fix the dictionary `L <-> q` (and `L = 1 <-> q = 0` in either convention).
*Source.* CP s.0.1 (i), s.1.6 (last paragraph); CP4 s.4 (B1)-(B10); RA s.1.3 (last item, the stronger remark).
*Status.* PROVED, R2 (RF s.1 CONFIRMED; RA s.1.1, s.1.3). *Dependencies.* Universality (F1) (= FSL Lemma U), (F2); Taylor data of D6;
`P_k in (1 - L^5)^-k Q[L]`; `-1/S_0''(L) = -L (1 - L^5)/5`.

**GAM-2** `prop:gam-forms` (prop). *Statement.* `c(z; L)` is also the normalised saddle series at `sigma = L` of (a) the integrand of
GAMMA s.2 (identical) and (b) POLY s.5's I-function summand `q^s Gamma(1 + 5s)/prod_(j=0..4) Gamma(1 + s - zeta^j/z)` in the flat
variable, all Gamma factors Stirling-expanded, with the logarithm of `q` chosen so that `L` is critical.
*Proof routes.* RF-C5 (direct, recommended): both flat phases have the same second derivative `5/sigma - Sum_j 1/(sigma - zeta^j) =
S_0''`, hence the same Taylor coefficients of order `>= 2`; the normalised amplitudes agree since `(zeta^j - sigma)^(1-2k) =
-(sigma - zeta^j)^(1-2k)`; `L` is critical for both by the choice of the logarithm. Alternative (CP s.1.7): the formal reflection
identity `log Gamma_St(y) + log Gamma_St(1 - y) = log 2pi + (y - 1/2)(log y - log(-y))` (PROVED in RL s.3.3 from FSL4 s.9.1 (iv) and the
oddness of `mu`), plus branch bookkeeping (`eps_j` odd, `e = eps + Sum eps_j` even, criticality forces `e = 0`; checked by RF s.1 and
RA s.1.1, not by RL).
*Source.* CP s.1.7; RF s.1 (s.1.7 item) and RF-C5; RA s.1.1; RL s.3.3. *Status.* PROVED, R2.
*Use.* Only to connect `c` with the hatted series `Phi^SP` (GAM-7). *Dependencies.* GAM-1; (L1a) for (b) in hatted variables.

## 2A.2 The local statement

**GAM-3** `lem:gam-deltasrc` (lem), "the delta lemma with a source" ((L2') in the sources). *Statement.* Variables `(s, m)` near
`0 in C^2`; `alpha0 in C^x`; phase `P(s) + alpha0 s m` with `P(s) = Sum_(k>=2) p_k s^k` holomorphic at 0 (EXACTLY this form: no pure
`m`-term and no term `s^k m` with `k >= 2`); amplitude `e^(cm) A(s, z)` with `c in C` and `A in O_0[[z]]`. Then `0` is a
non-degenerate critical point, `det H = -alpha0^2`, and, EXACTLY (no further corrections in `z`),
`N = exp(P(delta)/z) A(delta, z)`, `delta := -c z/alpha0`, both factors re-expanded in `z` (`P(delta)/z = Sum_(k>=2) p_k (-c/alpha0)^k z^(k-1)`).
For `c = 0` this is the delta lemma (`fsp` SP-L2 below) with the roles of multiplier and constrained variable exchanged: `N = A(0, z)`.
*Source.* RA s.3.2; LF s.2 (re-derivation, controls: a pure `m^3` term or an `s^2 m` term breaks it); LA s.1.4.
*Status.* PROVED, R2 (LF, LA); VERIFIED (RA2 s.4.5, `z^8`; LF2 L4, `z^12`). *Dependencies.* D5 only (Gaussian moments).
*Placement.* With Theorem LMI in the main text (one paragraph), or in the appendix; the writers choose.

**GAM-4** `thm:gam-lmi` (thm), "local Mellin inversion". *Statement.* Let `s*, alpha0, M0 in Z, g, F, a0, U*, lambda(L), S, A, c^F`
be as in D7. Then
* (i) each `c^F_n(L)`, `L in U*`, depends only on `L` (not on `lambda` or on the branch of `log a0`);
* (ii) each `c^F_n` extends holomorphically to `L = s*`;
* (iii) `c^F(z; s*) = exp(gt(delta)/z) F(s* + delta, z)/F(s*, 0)`, `delta = -M0 z/alpha0`, both factors re-expanded in `z`:
  `gt(delta)/z = Sum_(k>=2) g_k (-M0/alpha0)^k z^(k-1)`, `F(s* + delta, z) = Sum_n z^n Sum_(m+j=n) F_m^((j))(s*) (-M0/alpha0)^j/j!`.
In words: the saddle series, continued to the colour point, is the normalised value of the regular part of the integrand at the
first pole `sigma = s* + delta` of the singular factor `Gamma(M0 + a0/z)`.
*Scope* (LA-C6 = RS2 item 37; LA s.2.5): `g`, `F` holomorphic at `s*` and free of `lambda`; holomorphic dependence of `g`, `F` on `Q`
is allowed (the proof is verbatim; the value is that of the `Q = 0` data); `a0` exactly linear in `sigma`; `F(s*, 0) != 0` is a genuine
hypothesis (LA X9: simple poles otherwise). `M0 in C`: GAM-13.
*Source.* Statement RA s.3.1; proof RA s.3.3; (i) LF s.1; corrections LF-C1 (identification of the critical points written out),
LF-C2 (`Q^M0` drops out by `Nhat = N/N_0`, (F2) removes `lambda s*`), LF-C3 = LA-C3 (define `g_k`; the punctured disc `U*`), LA-C1
(`Omega_2` is a formal integrand), LA-C2 (full input list), LA-C6 (scope), LA-C7 (status); RS2 items 29-31, 34, 35, 37, 41.
*Status.* PROVED, R2 (LF "sound-with-corrections", LA "sound"); RS A8-A9; SYN5b K6. GEN s.3.3's "CONJECTURAL, SKETCH" is obsolete
(SYN5b O2; RS2 item 42: GEN's statement is RA's once `Lam(0)` is read as the constant term of `Lam(z)`).
*Dependencies* (LA-C2): (F1)-(F2); the holomorphic implicit function theorem (CP s.1.4; a formal IFT over `C[[Q]]` would also do,
LA s.2.1); (L1a) in its translation case and for `chi = exp`; (L1b) with the factor `beta`; (L1c)(i); (L3) at `a in C^x`, `M = M0`,
measure `dx` (FSL4 s.9.3 formal, or CP2 s.2.8 modulo DLMF); GAM-3; the `C[[z]]`-linearity of `N`.

## 2A.3 Theorem Gamma and its consequences

**GAM-5** `thm:gam-main` (thm), Theorem Gamma. *Statement.* For the series `c(z; L)` of D6:
* (i) `c_n in Q[L, L^-1, (1 - L^5)^-1]` (GAM-1);
* (ii) **(Gamma)**: no `c_n` has a pole at `L = 1` (`q = 0`), and `c(z; 1) = h(z)`, i.e. `c_n(1) = h_n` for every `n >= 0`;
* (iii) `c_n(zeta L) = zeta^-n c_n(L)`; hence no `c_n` has a pole at a fifth root of unity, `c_n in Q[L, L^-1]`, every exponent of `L`
  in `c_n` is `= 4n (mod 5)`, and `c(z; zeta^i) = h(zeta^-i z) = h_(w_i)(z)`;
* (iv) no `c_n` has a pole at `L = 0` (`q = oo`), `c(z; 0) = Ginf(z)` (GAM-9), and `c_n in Q[L]` (GAM-8).
*Source.* CP s.0.1 (CP's numbering: (ii) is the chart identity, (iii) is (Gamma), (iv) is Z/5); RA s.3.4; LF s.4.1-s.4.3; LA s.1.5;
RS2 s.6.1; M5 s.1. *Status.* (i)-(iii) PROVED by two independent proofs, each R2 (CP s.3: RF, RA; RA s.3: LF, LA); (iv) PROVED, R2
for `c(z; 0) = Ginf` (LF, LA) and R1 for `c_n in Q[L]` (LF s.4.3).
*Corrections.* Status K1 (GEN's "CONJECTURAL" superseded), O1, O18 (M4 s.1-s.2 superseded). CP's (ii) becomes a remark (GAM-11) by
M5 s.5 item 2. *Proposed proofs of record.* (ii) LMI at `s* = 1`, `alpha0 = -1`, `M0 = 0`; (iii) CP3 s.3.9 (or LMI at `s* = zeta^i`);
(iv) LMI at `s* = 0`, `alpha0 = 5`, `M0 = 1`. *Dependencies.* GAM-1, GAM-4, Lemma H (F3) and (L1a) for (iii), D3 (`T_k`), D6 (the
cancellation of the `sigma`-dependent powers of `z`, which identifies `c` with the summand's saddle series, RS2 s.6.1).
*At each of `L = 1, zeta^i, 0` exactly one Gamma factor degenerates* (`j = 0`, `j = i`, and the twist factor `Gamma(1 + 5s)`,
respectively), so LMI applies at all three points (LA s.1.5).

**GAM-6** `cor:gam-colours` (cor). *Statement.* `c^(i)(z; L) = c(z; zeta^i L) = c(zeta^-i z; L)`, i.e. `c^(i)_n(L) = zeta^(-in) c_n(L)`. The
statement (Gamma) at colour `i`, `c^(i)(z; 1) = h_(w_i)(z)`, is equivalent to (Gamma) at colour 0, because `w_i = zeta^i w_0` gives
`h_(w_i)(z) = h(zeta^-i z)`. So the covariance gives nothing beyond colour 0; on the gauge side it is `D_i(z) = D_0(zeta^-i z)` (UR s.4).
*Source.* CP s.0.2, CP3 s.3.9; RF2 s.4; RA s.1.3. *Status.* PROVED, R2. *Note.* Not needed by the chain: R0 (UR s.4) already gives
the colour relation (RF2 s.4; SYN5 S6).

**GAM-7** `prop:gam-hatflat` (prop), flat and hatted series. *Statement.* `c_n(L) = L^(4n) ph^SP_n(L^-5)` as rational functions of `L`;
hence `ph^SP_n(1) = c_n(1) = h_n` for every `n`, i.e. `Phi^SP(w; 1) = h(w)` (the identity called "(A2)" by POLY).
*Source.* CP3 s.3.10 with RF-C2; RL s.3.2 (Lemma H); SYN5 S3. *Status.* PROVED, R2 (RF2 s.5 item 2; RA s.1.3); VERIFIED (RF2 s.7
`rf_hat` (H3): `n <= 12` at `L = 2, -2/3, 3`).
*Correction.* RF-C2 (= SYN5b O5): the Jacobian `d sigma/dX = -L/X^2` of `sigma = L/X` is not constant; with the factor `X^2` of
`e^(S_1(L/X))` it gives POLY's `(1 - u X^5)^(-1/2)`, so it is part of the hatted amplitude; only its value at `X = 1` is removed by
normalisation. *Dependencies.* (L1a) for `sigma = L/X`; Lemma H with `mu = L^-4` (`z = L^-4 w`); (L1c)(i); GAM-2 (b); Theorem P(a)
(`ph^SP_n` is a polynomial) or GAM-5 (ii) for regularity at `u = 1` (RF2 s.5 item 2).
*Naming hazard* (SYN4 s.1.3): poly's "(A2)" is this identity; gamma's "(A2)" is an analytic step. Write neither name in the note.

**GAM-8** `cor:gam-poly` (cor). *Statement.* `c_n in Q[L]`: `c_n(L) = Sum_k a_(n,k) L^(4n-5k)` with every exponent `4n - 5k >= 0`.
With GAM-7: the top `u`-degree of `ph^SP_n` is `<= floor(4n/5)`, and `[u^(4m)] ph^SP_(5m) = c_(5m)(0)`. That `ph^SP_n` has no negative
power of `u` (equivalently `deg c_n <= 4n`) is Theorem P(a), or the regularity of the hatted series at `u = 0`; the chart does not give it.
*Source.* RA s.1.4 (R4), s.3.4; LF s.4.3; RS s.4; RS2 item 38. *Status.* PROVED, R1 (LF s.4.3); VERIFIED `n <= 30` (LF2 F).
*Corrections.* RF-C3 = SYN5b O6 (CP4 s.5 item 7 "PROVED elsewhere per SYN4 s.1.1" overclaims: polynomiality needs P(a), s.3.10 and the
degree bound, or the chart at `L = 0`); RS2 item 38 ("without Theorem P(a)" covers `c_n in Q[L]` only).
*Dependencies.* GAM-5 (i), (iii), (iv) at `L = 0`; GAM-7 for the `ph^SP` form.

**GAM-9** `thm:gam-ginf` (thm or prop), the orbifold point. *Statement.* `c(z; 0) = Ginf(z)`. Equivalently: at `s* = 0`, `alpha0 = 5`,
`M0 = 1` (`a0 = 5 sigma`, `delta = -z/5`, `Q = (-q)^(-1/5)`), the LMI value is the normalised value of `prod_(j=0..4) Gamma_St((zeta^j - sigma)/z)`
at `sigma = -z/5`, i.e. of `prod_j Gamma_St(zeta^j/z + 1/5)`, normalised by division by `(2 pi)^(5/2) z^(3/2) e^(g(0)/z) e^(-g'(0)/5) F(0, 0)`.
*Source.* RA s.3.4 (orbifold bullet); LF s.4.2 with LF-C5; LA s.1.5, s.2.3, LA-C5; RS A7; SYN5b K7.
*Status.* PROVED, R2 (LF, LA), purely formal with the formal DLMF 5.11.8 (RS A7; `fsp` SP-5118); VERIFIED (LF2 Q4 chart at `L = 0`
through `z^60`, D3 value step `z^80`; RA2 O1 `z^18`; GEN2 s.5.1 A `n <= 20`; LA2 Q4).
*Corrections.* LF-C5 (state the normalisation, as above); LA-C5 (the formal lemma replaces DLMF 5.11.8 via uniqueness of asymptotic
expansions, which is also valid, LF s.4.2); SYN5b O14 (RA's "given the inputs of Theorem D" asks for more than is needed).
*Dependencies.* GAM-4 at `s* = 0`; SP-5118; FSL4 s.9.1 (iv); `B_k(theta + 1) - B_k(theta) = k theta^(k-1)` (DLMF 24.4.1, itself from the
generating function 24.2.3, LA s.2.3).

**GAM-10** `cor:gam-free` (cor), the free coefficients. *Statement.* `f^SP_m = [u^(4m)] ph^SP_(5m) = c_(5m)(0) = [z^(5m)] Ginf` for all
`m >= 1`. Given Theorem P(b) (`Phi^SP = Phi^can`): `Sum_m f^can_m z^(5m) = Ginf(z)`; no (QRR) and no `D = 1` is needed for this. The
equality `f^QRR = f^can` needs `D = 1`.
*Source.* RA s.3.4; LF s.4.3; SYN5 s.3 (last item); SYN5b K7; RS s.4; M5 s.1. *Status.* PROVED; the step to `f^can` R1 (LF s.4.3).
The closed form was CONJECTURAL in SYN4 s.1.1 and M4 s.3 (VERIFIED `m <= 100` at two primes, `m <= 25` exactly).
*Prior art.* GJR arXiv:1812.11908v1 Lemma 7.7, eq. (37) (`m = 5`, `k = 1`, `lambda = 1`) proves the same formula for the normalised SP
series of their LG integral at the orbifold point (LIT2 s.3.6). The identification of GJR's series with `c(z; L)` is written nowhere
(LIT2 s.3.6); the note proves `c(z; 0) = Ginf` directly and cites GJR as prior art for the formula (SYN5b K7: the transfer route of
LIT3 s.4 item 2 is superseded). *Dependencies.* GAM-8, GAM-9; Theorem P(b) (other digest).

## 2A.4 Remarks (second proof, by-products, sharpness)

**GAM-11** `rmk:gam-sevenvar` (rmk), the 7-variable chart. *Statement.* For every `L in Lam`: `c(z; L) = Nhat_7(z; L) = Nch(z; L)` (D8);
`Nch_n in Q(zeta)[L, L^-1, (zeta^j - L)^-1 : j = 1..4]`, regular at `L = 1`; `Nch(z; 1) = exp(Sum_k b_k T_k z^(2k-1)) = h(z)`.
Data: `det H_7 = 1/(L^2 (1 - L^5)^2)`; order (a) inner block `diag(-1/(5L), -1/(zeta^j - L))`, Schur complement `S_0''(L)`; order (b)
block `[[0, -1/(1-L)], [-1/(1-L), -1/(1-L)]]`, `det = -1/(1-L)^2`, Schur complement `= Hess psi(p_ch(L))` entrywise;
`det H_ch = -1/(L^2 Lam1(L)^2)`; `beta = X_0(w)/X_0(w_*)` ("the Jacobian `X_0` cancels the `1/x_0` of `dx_0/x_0`").
*Correction RF-C1* (= SYN5b O4): CP3 s.3.2's "`det H_ch` ... vanishes only at `L = 0, zeta^j`" is false; read "`1/det H_ch = -L^2 Lam1(L)^2`
vanishes only at `L = 0, zeta, .., zeta^4`, the `L` at which `p_ch(L)` leaves the chart; for every other `L`, `H_ch` is invertible, and
`det H_ch(1) = -1/25`". *Source.* CP s.0.1 (ii); CP3 s.3.1-s.3.6; RF s.3; RA s.1.3; FSL3 s.8. *Status.* PROVED, R2.
*Warning* (FSL3 s.8 warning 1; SYN5b K10, no conflict): the 7-variable SP is not defined at `L = 1`; both orders are used on `Lam` only,
and the value is read off the chart, which is defined at `L = 1`. *Placement.* Remark after Theorem Gamma (M5 s.5 item 2).
*Dependencies.* (L0), (L1a), (L1b) twice, (L1c)(i), (L2), (L3) with `M = 1` (`dy`) and `M = 0` (`dx_j/x_j`), (F1)-(F3).

**GAM-12** `rmk:gam-Iexp` (rmk), coefficientwise agreement with the I-function. *Statement.* In `Q((z))[[q]]`,
`L e^(Dpsi(q)/z) c(z; L(q)) = h(z) Sum_(d>=0) q^d I_d(z)`, `L(q) = (1 - 5^5 q)^(-1/5)` (the series with `L(0) = 1`),
`Dpsi(q) := psi_q(p_ch(q)) - psi_0(p_ch(1))`, `D_q Dpsi = L - 1`; so `log q + Dpsi` is, up to a constant, the colour-0 canonical
coordinate (`D_q u_0 = L`). *Source.* CP3 s.3.7, CP4 (B10); RF s.3.7; RA s.1.3. *Status.* PROVED via (L1c)(iii), R2; not used.
*Open.* CP3 s.3.7 calls it "the form (Gamma3)" but says the identification with arXiv:0710.0889 eq. (5) was not re-checked; no source
writes it out (needs a `\draftnote` if mentioned). *Dependencies.* (L1c)(iii) (Corollary T'), (L0), (L3) with `M = 5d + 1` and `M = -d`.

**GAM-13** `rmk:gam-Mcomplex` (rmk). *Statement.* LMI holds for every `M0 in C` with `A := a0^(M0 - 1/2) exp(Sum_n beta_n(M0) (z/a0)^n) F`
(RA's proof verbatim, with SP-L3C in step 1; `Q^M0 := e^(-M0 lambda/alpha0)` a constant; GAM-3 allows any `c`). The naive analogue
(elementary part of `Gamma(M0 + w)` with the Bernoulli part of `Gamma(w)`) fails: pole at `z^1`, coefficient `(beta_1(0) - beta_1(M0))/alpha0`.
*Source.* LA s.2.3 (Corollary); LF-C6; RS A3. *Status.* PROVED by LA; checked by hand and VERIFIED only by the RS critic (RS2 s.6.3;
`rc_l3c`); not used by the chain. *Conflict.* LF-C6 (SKETCH, for RA2 s.5's sentence) vs LA (PROVED) -> RS A3: LA; RS2 item 39.

**GAM-14** `rmk:gam-blind` (rmk). *Statement.* (Gamma) is blind to the measure of the twist factor: with `Gamma(5s)` (measure `dy/y`)
in place of `Gamma(1 + 5s)` (`dy`) the value at `L = 1` is still `h` (`Poch_0 = Poch_1 = 1`; in LMI the extra factor `5 sigma` is `z`-free
at `sigma = 1`), while `c_n(L)` changes at `L != 1` exactly for `n >= 2`. The correct `dy` enters the chain through Theorem P(b) and GAM-7.
At `L = 0` the twist factor is the singular one: `M0 = 1` gives `Ginf`, `M0 = 0` gives `exp(5 Sum_(5 | 2k-1) b_k z^(2k-1)) != Ginf`.
*Source.* RF s.1 (s.1.7 remark), RF2 s.6 (remark); LA s.2.4; LA2 s.6 (write-up note); RS A10; RS2 item 40. *Status.* PROVED (from (L3)
and LMI); VERIFIED (RF2 `rf_sigma`, `z^6`; LA2 Q2, Q5). *Use.* A sentence in the note, so that (Gamma) is not read as a test of `dy`.
