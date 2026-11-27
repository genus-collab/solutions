# DIGEST_D2, part 3: s.2 Statements, end of part A (remarks); part B (formal stationary phase, for the appendix)

Continuation of `DIGEST_D2_2.md`. (QRR) statements: `DIGEST_D2_4.md` s.2C.

## 2A.5 Remarks: sharpness, generality, prior art, interface

**GAM-15** `rmk:gam-necessity` (rmk), what the hypotheses of LMI do. *Supported statement* (RA s.3.4 "What this says"; M5 s.3.5;
GEN2 s.5 "Reading"): regularity at `q = 0` is a property of the degenerating factor alone; it must be an exact `Gamma_St` form and the
only factor degenerating there; the regular factors and the twist enter only the value. The "if" direction is LMI. The "only if"
direction is EVIDENCE, except (a):
* (a) PROVED (GEN s.3.2 Lemma 4; R0, GEN was not refereed; read by the SYN5 critic): if the argument of a twist factor
  `Gamma(e + l s + d/z)` vanishes at `sigma = lambda_i` (`l lambda_i + d = 0`, `l >= 1`, `lambda_j != lambda_i`), no critical point tends to
  `lambda_i` as `q -> 0`: no colour exists.
* (b) The one-Gamma identity (LF s.1, "a structural check"; derivation, not labelled as a lemma): for `g = 0`, `F = 1`, (L1a) with
  `sigma - s* = (L - s*) t` and (F3) with `mu = a0(L)` give `c^F(z; L) = C(z/a0(L))` for one universal series `C` (depending on `M0`);
  regularity at `s*` forces `C = 1`, which is LMI's value. "The formal shadow of Mellin inversion `(1/2 pi i) Int Gamma(M0 + a) x^(-a) da =
  x^M0 e^(-x)`." By hand `N_1 = 1/12 - 3/8 + 1/2 - 5/24 = 0` (`M0 = 0`). VERIFIED `C = 1` through `z^12`, `M0 = -2..3` (LF2 L2); GEN s.1
  "layers": the leading pole coefficients `ct_n(0)` come from the singular factor alone and vanish (`n <= 20`).
* (c) VERIFIED only (evidence appendix; GEN2 s.5.3, RA2 s.4.4, LA s.2.2-s.2.5, s.3, LA2 s.4.3): two degenerate factors (coincident
  weights; `Gamma(a0/z)^2`, `Gamma(a0/z) Gamma(-a0/z)`, `Gamma(2 a0/z)`, `Gamma(1 + a0/z)`; two slopes at `q = oo`) give poles, of order
  exactly `n` at `z^n` for the square (leading coefficients `((2n-1)!!)^2/(n! 16^n) (-1/alpha0)^n`, Hankel's `K_0` coefficients "used only
  as an explanation"), order `2n` for `Gamma(a0/z) Gamma(-a0/z)`; a perturbed singular `b_k` first gives a pole at `z^(2k-1)` of order
  `2k - 1` (coefficient `eta/alpha0^(n0)` for a perturbation `eta` of the `(z/a0)^(n0)` coefficient, X6); a wrong exponent `a0^E`,
  `E != M0 - 1/2`: pole at `z^1`, coefficient `(1/24 + beta_1(M0) - E^2/2)/alpha0` ("derived from the five order-z Wick terms", LA2 s.4.3);
  `F = 1/(sigma - s*)`: poles of order `n`, coefficient `n! (-1/alpha0)^n` (observed; `E_1`, DLMF 6.12.1, an explanation, not a proof);
  `F(s*, 0) = 0`, `M0 != 0`: simple poles; a non-linear singular argument: dependence on `log Q`, non-holomorphic limits (X4); a
  `lambda`-dependent `F`: logarithmic divergence (X5). Perturbing a REGULAR factor never creates a pole and changes the value exactly as
  LMI predicts (RA2 I1-I4; GEN s.4.2; LA2 P). Two "falsifiers" of RA's first draft were covered cases (`+1/2` exponent = `Gamma(1 + w) =
  w Gamma(w)`; `1/Gamma` = reflection), RA s.2, RA2 s.7.
*Placement.* A remark after LMI stating (a) and (b) with proofs, and pointing to the evidence appendix for (c) (rule 3: no proof cites it).

**GAM-16** `rmk:gam-general` (rmk), generality. LMI applies verbatim to any Mellin-Barnes summand with exactly one degenerate factor
with a linear argument and `lambda`-free regular data (LA2 s.5: GEN's T2, T3 and GEN2 s.5.1 A-E are instances of LMI). Cases computed
(GEN s.0, s.4, GEN2 s.5.1-s.5.2; VERIFIED): O(5)-twisted `P^4` with four sets of distinct rational weights, all 20 colours, value
`h_(w_i)` with `w_i = {lambda_i - lambda_j} U {-5 lambda_i}` (`n <= 24`); `P^1`, `P^2`, local `P^2` (two lifts), `O(-1) + O(-1)` over `P^1`
(two lifts), `O(3)`, `O(4)`, `O(6)` over `P^4` and `O(4)` with a shifted lift: 12 theories, 40 colours (`n <= 20`); the `q = oo` charts
(twisted `P^4`: `exp(Sum_n beta_n(1/5) z^n Sum_j lambda_j^-n)`; local `P^2`: no shift); the shift rule `h_f = e_f - c_f e_*/c_*` (a
corollary of LMI and SP-5118, LA2 s.5); an orbifold-type factor with `c_* = -2`; random regular WKB data. The QRR sign rule seen in
these cases (GEN s.3.1, s.4.1): an Euler-class twist contributes minus its fibre weight, an inverse Euler class plus its fibre weight;
"I did not re-derive it from the QRR papers" (GEN s.4.1). Hypotheses NOT used (GEN s.3.2): the CY condition, `Sum lambda_j = 0`, `|c_*| = 1`,
`e_* = 0`, the sign of `q`, the Gamma form of the regular factors. *Status.* The summand/I-function link is GEN Lemma 2 (PROVED as a
residue identity, R0; I-function formulas quoted without statement numbers, SYN5b O3); the regularity is LMI; the QRR of those
theories is not quoted anywhere. *Placement.* One remark; the cases in the evidence appendix. GEN's labels S2-S6 (CONJECTURAL) predate LMI.

**GAM-17** `rmk:gam-prior` (rmk), prior art (details and statement numbers in `DIGEST_D2_6.md` s.6). Givental arXiv:math/0108100v2 s.10
gives the untwisted mechanism at `Q = 0` for equivariant `CP^n`: the chart `x_0 = Q/(x_1..x_n)`, factorisation into one-variable integrals,
Lemma 10.4; continuity in `Q` is asserted (LIT2 s.3.1; RF2 s.8). Differences (CP3 s.3.8): the twist (variable `y`, factor `Gamma(1 + 5/z)`,
weight `-5`); the object of (Gamma) is the one-variable series `c(z; L)`; the constraint comes from the Mellin-Barnes variable `sigma`
through the delta lemma (`sigma` is Givental's Lagrange multiplier `P`, p. 14); the continuation acts on coefficients. GJR
arXiv:1812.11908v1 Remark 7.2 states the analogue of (Gamma3) for their LG integral without proof; their Lemma 7.7, eq. (37) proves the
`Ginf` formula for it. Iritani arXiv:0903.1463v3 Lemma 4.20: an exact `q`-expansion at `z = 1` in the analogous chart. "No paper read
proves (Gamma). The write-up should cite these and claim no more novelty" (M5 s.3.6; LIT3 s.4 item 3; SYN5b s.8 item 2).

**GAM-18** (no label; interface with the chain). Theorem Gamma supplies `c_n(1) = h_n`; GAM-7 turns it into `Phi^SP(w; 1) = h(w)`; with
Theorem P(b) and the glue `D_0(z) = c(z; 1)/h(z)` (SYN5 S4-S6; other digests) this gives `D_0 = 1`. CP s.0.3/s.3.11 ("Theorem D") is a
valid but longer deduction whose hypothesis list is incomplete and contains an unnecessary item: it omits (CYC), (S2), the relabelling,
and includes R2 = Lemma N (SYN5b O8, with K3, K4, K5, K12; RF's "complete" in RF s.0 is likewise corrected). Use SYN5 s.1-s.2 with RS2
s.6.2. Not transcribed here.

## 2B. Formal stationary phase (proposed appendix, prefix `fsp`)

**Two complete written sets exist**: CP s.1-s.2 (Wick functional `N`, lemmas (S1), (PG), (L0), (L1a)-(L1c), (L2), (L3), facts
(F1)-(F4)) and FSL s.1-s.10 (weighted ring, `E_zC`, Lemma U, (G1)-(G5'), FL1-FL4, IFT, Lemma T, Corollary T', L1a with the Piola
identity, L1b, L2, L3 with a formal Proof A, the formal residue). Both refereed (RF s.2; RL s.1). SYN5b s.8 item 2 (vii): the appendix is
"FSL, with RL's Lemma H". Dictionary (RL s.3.1 table; SYN5b K11):

| CP name | FSL / RL name | proposed label | users in the chain (and remarks) |
|---|---|---|---|
| (F1) | Lemma U (FSL "(U)": rename, clash with Theorem U) | `lem:fsp-univ` | GAM-1, GAM-4 (LA-C2), (L1c)(i), L3-C |
| (F2) | implicit in FSL Def. 1.2 | `lem:fsp-const` | GAM-4 step 2 |
| (F3) | RL Lemma H (FSL4 (H) is a special case) | `lem:fsp-homog` | GAM-5 (iii), GAM-7, L3 step (1), P(b) step 4, P(a) |
| (S1) | FL1 | `lem:fsp-ibp` | (L1a) proof, Lemma T, L3 A1, P(b) step 2, P(a) |
| (PG) | (G5), (G5') | `lem:fsp-pg` | (L1b), Corollary T', FL4 |
| (L0) | (G4) | `lem:fsp-prod` | (L1b); GAM-11, GAM-12 |
| CP s.1.4 (hol. IFT) | formal IFT (FSL2 s.3.2) | `lem:fsp-ift` | GAM-4 step 2, (L1b), (L2), T' |
| - | Lemma T (a), (b) | `lem:fsp-transl` | Corollary T'; P(a) shift `sigma -> sigma + z` (T(b)) |
| (L1c)(iii) | Corollary T' = L1c(iii) | `cor:fsp-movebase` | GAM-12; L1c(iv) |
| (L1c)(i), (ii) | L1c (i), (ii) | `lem:fsp-param` | GAM-1, GAM-4, GAM-7, GAM-11 |
| - | L1c (iv) (T-z) | in `lem:fsp-param` | L3 Proof A step A2; gamma's (S4); P(a) |
| (L1c)(iii') | FL4 | `lem:fsp-dq` | P(a) (HPF); PB's (S2) (PS-C1) |
| (L1a) | L1a (+ Piola identity) | `lem:fsp-cov`, `lem:fsp-piola` | GAM-2, GAM-4, GAM-5 (iii), GAM-7, (L1b), L3 |
| - | FL3 for the SP (FSL2 s.5.4) | `lem:fsp-cfu` | P(b) step 3 |
| (L1b) (Fub) | L1b | `lem:fsp-fubini` | GAM-4 step 1; GAM-11 |
| (L2) | L2 | `lem:fsp-delta` | GAM-11 only (LMI uses GAM-3) |
| (L3) | L3 (i)-(v) | `lem:fsp-stirling` | GAM-4 step 1; GAM-11, GAM-12 |
| (functional equation) | L3 (iv) | `lem:fsp-stirfe` | L3 Proof A; SP-5118; reflection; P(a) |
| - | formal DLMF 5.11.8 (LF s.4.2, LA s.2.3) | `lem:fsp-dlmf5118` | GAM-9 |
| - | Lemma L3-C (LA s.2.3) | `lem:fsp-L3C` | GAM-13 only |
| reflection (CP s.1.7) | RL s.3.3 | `lem:fsp-refl` | GAM-2 (alternative route) |
| residue items of (FL) | Lemma 10.1, FL2, FL3 (FSL4 s.10) | `lem:fsp-res` | P(b) steps 3, 6-7 |
| (F4) | FSL s.1.4 "Meaning" | `rmk:fsp-analytic` | none (remark; never in a proof) |

*Numbering hazard* (RL s.3.4; SYN5b K11): PB_AUDIT s.2.0 numbers its own FL1-FL7; PB FL2 (rank-one quotient) and PB FL3 (parameter
derivatives) are NOT FSL's FL2, FL3. Dictionary: PB FL1 = FSL FL1; PB FL2 = FSL L1a and FL3-SP; PB FL3 = FSL L1c(iii)/FL4; PB FL4-FL5 =
FSL FL2-FL3; PB FL6 = FSL L3(iv); PB FL7 = FSL Lemma T(b). The note should use descriptive names only.
*Ledger* (RL-C6 = SYN5b O16, O10, K11): the formal input of Theorem P is not only the poly referee's four lemmas (exact forms, changes of
variables termwise in `u`, `D_q`); it also uses L3 (iv), Lemma T(b), Lemma H and the relabelling `V = s Res_X(..)` (FSL4 s.10.4; PB s.2.7).

**Statements (final form; proofs in `DIGEST_D2_5.md`)**

* **SP-DEF** `def:fsp-wick` (defn). D5; FSL s.1.2 (weight `wt(t^alpha z^j) = |alpha|/2 + j`, ring `R_n(A)`, Lemma 1.1: a ring, `d/dt_i`
  continuous, `exp(P)` converges if every term of `P` has weight `>= 1/2`); FSL s.1.3 (`<t^alpha>_C := d_xi^alpha exp(xi^T C xi/2)|_0`, the Wick
  recursion (W), `E_zC[t^alpha z^j] := <t^alpha>_C z^(j + |alpha|/2)`); FSL Def. 1.2 `SP^norm_c[F; a] := E_zC[exp(F_3/z) a(c + t)]` (= CP's
  `N`, RL-C7; the note's normalised series `Nhat` divides by its constant term `a(c)`); Def. 1.3 full SP (symbols). Over any commutative
  Q-algebra `A` with `H in GL_n(A)`. PROVED (Lemma 1.1), R1 (RL).
* **SP-U** `lem:fsp-univ` (lem). `[z^k] SP^norm` is a polynomial with rational coefficients in finitely many Taylor coefficients of the
  phase, of the amplitude and in the entries of `C`; for an amplitude `Sum_m A_m z^m`, the orders are `3..2k+2` for the phase and
  `<= 2(k - m)` for `A_m` (CP (F1), sharp; FSL's `3..6(k - j0)` is a correct, non-sharp bound, RL s.1). (U1) ring maps commute with
  `SP^norm`; (U2) holomorphic or formal dependence on parameters is inherited. CP s.1.2 (F1); FSL s.1.5. PROVED, R2 (RF s.1, RL s.1).
* **SP-CONST** (F2): adding a constant to the phase changes only the critical value. Trivial.
* **SP-H** `lem:fsp-homog` (lem), homogeneity. For `mu` a unit: `N[mu Phi, a(., z)](z) = N[Phi, a(., mu z')](z')` at `z' = z/mu`; symbols
  `[e^(mu Phi(c)/z)] (2 pi z)^(n/2) [det(-mu H)^(-1/2)] = [e^(Phi(c)/z')] (2 pi z')^(n/2) [det(-H)^(-1/2)]` with `det(-mu H)^(1/2) := mu^(n/2)
  det(-H)^(1/2)`. CP s.1.2 (F3) (CONFIRMED RF s.1, RA s.1.1); RL s.3.2 Lemma H (missing statement in FSL, RL-C1). PROVED.
* **SP-IBP** `lem:fsp-ibp` (lem). `N[Phi, Sum_i ((d_i Phi) h_i + z d_i h_i); p] = 0` for amplitudes `h_i` (FSL: any `h_i` with
  `e^(F_3/z) h_i in R_n(A)`, and every sum of such expressions that converges in `R_n(A)` after multiplication by `e^(F_3/z)`). CP2 s.2.1;
  FSL2 s.3.1 (with the Remark on the series `g(sigma) - g(sigma + z)` of Theorem P(a): weight `>= (k+1)/2`). PROVED, R2.
* **SP-PG** `lem:fsp-pg` (lem). CP form: `R` complete for an ideal `I`, `H` invertible symmetric, `delta` symmetric with entries in `I`,
  `C_delta := -(H + delta)^-1`: `<exp((1/2) eta^T delta eta) P>_C = det(1 + H^-1 delta)^(-1/2) <P>_(C_delta)`, `I`-adically, the root by the
  binomial series; infinitesimal form `d/dt [det(-H_t)^(-1/2) <P>_(C_t)] = det(-H_t)^(-1/2) <(1/2) eta^T (dH_t/dt) eta P>_(C_t)`. FSL form
  (G5'): `E_(zC_0)[e^(t^T K_1 t/(2z)) P] = (det K/det H_0)^(-1/2) E_(-z K^-1)[P]`, `K = H_0 + K_1`, `K_1` topologically nilpotent. RL-C4
  (wording of the "ideal"), RL-C5 (state the case `K_1 = K(y) - H_xx` depending on outer variables, `P in R_(n+m)(A)`, used in L1b;
  same proof). CP2 s.2.2; FSL s.2. PROVED, R2.
* **SP-PROD** `lem:fsp-prod` (lem). Separated phases and amplitudes: `N` factorises (block-diagonal covariance). CP2 s.2.3; FSL (G4). R2.
* **SP-IFT** `lem:fsp-ift` (lem). `B` complete and separated for `I`, `F in B[[t]]`, `grad F(0) in I^n`, `Hess F(0) in GL_n(B)`: a unique
  `c in I^n` with `grad F(c) = 0`; `Hess F(c)` invertible. RL-C4: "the bracket is a matrix with entries in `I`, times `(c - c')`".
  FSL2 s.3.2. PROVED, R1. The holomorphic IFT (CP s.1.4) is standard.
* **SP-T** `lem:fsp-transl` (lem), Lemma T. `p` small: (a) `p in (qA)^n` over `A = A_0[[q]]`, or (b) `p in (zA[[z]])^n`. With
  `b_p(t) := exp((F(c+t+p) - F(c+t))/z) a(c+t+p)`: `SP^norm_c[F; b_p] = SP^norm_c[F; a]` (`b_p` admissible). FSL2 s.3.2. PROVED, R1.
* **SP-T'** `cor:fsp-movebase` (cor), Corollary T'. `F = F_0 + delta`, `delta in q A_0[[t]][[q]]`: unique critical point `c in c_0 + (qA)^n`,
  and `SP^norm_(c_0)[F_0; e^(delta/z) a] = e^((F(c) - F_0(c_0))/z) (Delta(c)/Delta_0)^(-1/2) SP^norm_c[F; a]` in `A_0((z))[[q]]`,
  `Delta = det(-Hess)`. = CP (L1c)(iii) (iterated derivative form). FSL2 s.3.2; CP2 s.2.6 (iii). PROVED, R2.
* **SP-PAR** `lem:fsp-param` (lem), L1c. (i) holomorphic (or formal) data with a non-degenerate critical point at `lam_0`: the critical
  points form a unique holomorphic family, and every coefficient of `SP^norm`, the critical value and `Delta` are holomorphic; rational
  data give rational coefficients with poles only where the Taylor data have poles or `det H = 0` (CP form); (ii) Taylor expansion and
  evaluation at `lam_0` commute with `SP^norm`; (iii) = SP-T'; (iv) T-z: for `g` with `g(c) = 0`, the series
  `S(eps; z) := e^((F_eps(c(eps)) - F(c))/z) (Delta_eps/Delta)^(-1/2) SP^norm_(c(eps))[F + eps g; a]` lies in the weighted ring with `eps` of
  weight 1 (every term of `eps g/z` has weight `>= 1/2`, RL-C4), and `[S(eps; z)]_(eps = z) = SP^norm_c[F; e^g a]`. CP2 s.2.6; FSL2 s.4.
  PROVED: (i)-(iii) R2; (iv) R1 (RL s.2: "the one genuinely new device"; VERIFIED in its L3 instance, RL2 s.7.2).
* **SP-DQ** `lem:fsp-dq` (lem), FL4. `dS/dq = e^((..)/z) (Delta(q)/Delta(q_0))^(-1/2) SP^norm_(c(q))[F; z^-1 (d_q F) a + d_q a]`; in full-SP
  language `d/dq SP[e^(F/z) a dx] = SP[d_q(e^(F/z) a) dx]`; the iterated clause `P(z D_q) SP[..] = SP[P(sigma) ..]` (when `D_q F = sigma`)
  needs `D_q a = 0` (RL-C2; true for the Mellin-Barnes amplitude). FSL2 s.4.5; CP2 s.2.6 (iii'). PROVED, R2.
* **SP-COV** `lem:fsp-cov` (lem), L1a. `phi(c') = c`, `J := D phi(c')` invertible (analytic, formal or parameter-dependent maps):
  `SP^norm_c[F; a] = SP^norm_(c')[F o phi; (a o phi) det D phi/det J]`, Hessian `J^T H J`; full form with `det(-J^T H J)^(1/2) := det(J)
  det(-H)^(1/2)`. CP form: `N[Phi o chi, (A o chi) J; p'] = J(p') N[Phi, A; p]`; `Nhat` invariant. RL-C3 (`P_s := G o psi_s - G(0) - y^T H y/2`).
  CP2 s.2.4; FSL2 s.5 (+ Piola identity `Sum_k d_k adj(D psi)_(ki) = 0`, `lem:fsp-piola`). PROVED, R2.
* **SP-CFU** `lem:fsp-cfu` (lem). L1a over `Q[[u]]` for the CFU map `Y = X psi(X^5; u)`, `psi in 1 + u Q[X^5][[u]]`, at `X = 1` (and
  `X = zeta^-i` over `Q(zeta)[[u]]`). FSL2 s.5.4. PROVED, R1. User: P(b) step 3 (other digest).
* **SP-FUB** `lem:fsp-fubini` (lem), L1b. `x = (v, w)`, `H_vv` invertible at `p`: (i) unique inner critical point `v_c(w)`; outer phase
  `phi(w) = Phi(v_c(w), w)`, `Hess phi(w_0) = S` (Schur complement), `det H = det H_vv det S`; (ii) the inner series `B(w, z)` has
  holomorphic (formal) coefficients; (iii) `N[Phi, A; p] = N[phi, beta B; w_0]`, `beta(w) := (det H_vv(w)/det H_vv(w_0))^(-1/2)`, `beta(w_0) = 1`;
  symbols `det(-H)^(1/2) = det(-H_vv)^(1/2) det(-S)^(1/2)`. Two splittings compute the same total SP. RL-C3 (`G_3 := G - G(0) - y^T H_yy y/2`),
  RL-C5. CP2 s.2.5; FSL3 s.6. PROVED, R2.
* **SP-L2** `lem:fsp-delta` (lem). `Phi = sigma G(x_0, w) + f(x_0, w)` EXACTLY linear in `sigma`, amplitude free of `sigma`,
  `G(X_*, w_*) = 0`, `d_(x_0) G(X_*, w_*) != 0`: `(sigma_c(w), X_0(w))` is a non-degenerate critical point, `det = -(d_(x_0) G)^2`, critical
  value `f(X_0(w), w)`, and `N = A(X_0(w), w, z)` EXACTLY; with SP-FUB the outer amplitude is `A(X_0(w), w, z) g_c/G_t(X_0(w), w)` ("restriction
  to `G = 0` with the Jacobian `1/G_t`, to all orders"). CP2 s.2.7; FSL3 s.7. PROVED, R2 (RF, RA for CP's proof; RL for FSL's).
* **SP-FE** `lem:fsp-stirfe` (lem). `rho(Y) := (Y + 1/2) log(1 + 1/Y) - 1 + mu(Y+1) - mu(Y) = 0` in `Q[[1/Y]]`, `mu(Y) := Sum_k b_k Y^(1-2k)`,
  i.e. `Gamma_St(Y + 1) = Y Gamma_St(Y)` formally; equivalently `St_(m+1)(w) = (1 + m w) St_m(w)`, `m in Z`. FSL4 s.9.1 (iv), s.9.2. PROVED, R1.
* **SP-L3** `lem:fsp-stirling` (lem). `a in A^x` (a unit), `m in Z`: (i) `SP^norm_a[a log x - x; (x/a)^(m-1)] = St_m(z/a)`, equivalently
  (CP form) `N[a log x - x, x^(M-1); a] = a^(M-1) Poch_M(z/a) exp(Sum_k b_k (z/a)^(2k-1))`; (ii) symbol form
  `SP_a[x^(a/z + m) e^(-x/z) dx/x] = [a^(a/z) e^(-a/z)] (2 pi z)^(1/2) a^(m - 1/2) St_m(z/a)`; measure `dx` = `m + 1`; (iii) this is the expansion
  of `z^(a/z + m) Gamma_St(a/z + m)`; (iv) SP-FE; (v) `SP_a[x^(a/z+m+1) e^(-x/z) dx/x] = (a + m z) SP_a[x^(a/z+m) e^(-x/z) dx/x]`. Insertion
  `x^m` with `dx/x` is `M = m`, with `dx` it is `M = m + 1`. CP2 s.2.8; FSL4 s.9.1-s.9.4. PROVED, R2 (formal Proof A: FSL4 s.9.3, RL).
* **SP-5118** `lem:fsp-dlmf5118` (lem). In `Q[theta][[1/w]]`: `(w + theta - 1/2) log(1 + theta/w) - theta + Sum_k b_k (w + theta)^(1-2k) =
  Sum_(n>=1) beta_n(theta) w^-n`, i.e. `log Gamma_St(w + theta) - [(w + theta - 1/2) log w - w + (1/2) log 2pi] = Sum_n beta_n(theta) w^-n`.
  LF s.4.2; LA s.2.3 (Lemma F5.11.8). PROVED twice independently (same argument), checked by hand and VERIFIED by the RS critic (RS A7,
  RS2 s.6.3). Formal replacement of DLMF 5.11.8.
* **SP-L3C** `lem:fsp-L3C` (lem). For `M in C`, `a in A^x`: `N[a log x - x, x^(M-1); a] = a^(M-1) exp(Sum_n beta_n(M) (z/a)^n)`; for `M in Z`
  the CP form. LA s.2.3. PROVED (LA), checked by the RS critic only (RS A3). Not used by the chain.
* **SP-REFL** `lem:fsp-refl` (lem). `log Gamma_St(y) + log Gamma_St(1 - y) = log 2pi + (y - 1/2)(log y - log(-y))` as formal Stirling
  expressions (the `mu`-terms cancel, `1 - 2k` odd). RL s.3.3 (from SP-FE at `Y = -y`). PROVED (RL). Branch part not a formal lemma.
* **SP-RES** `lem:fsp-res` (lem). The ring `L` of FSL4 s.10.1 (Lemma 10.1), `Res: L -> K((v))`; FL2 (`Res` kills `d/dX`-exact elements);
  FL3 (invariance under `X -> X + eps`, `eps in u K[X, 1/X][[u]]`, termwise in `u`). FSL4 s.10.1-s.10.3. PROVED, R1 ("trivial in `L`", RL).
  NOT covered: the relabelling `V = s Res_X(..)` (FSL4 s.10.4 SKETCH; PB s.2.7, PROVED R1 by PS) and gamma's (R3) (FSL5 s.12; not needed
  by P(b)). User: P(b) (other digest).
* **SP-ANAL** `rmk:fsp-analytic` (rmk). For real data with negative definite `H`, the full SP is the Laplace asymptotic expansion over a
  small ball. Used in no proof. CP (F4) cites "Hoermander, ALPDO I, Thm 7.7.5" FROM MEMORY, not re-read (CP s.1.2; LIT3 s.5): do not pin
  it; FSL s.1.4 states the meaning without a citation.
