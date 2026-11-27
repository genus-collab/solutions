# Digest D1: the connection problem (normalisation, polynomial structure, Theorem P, the glue)

Design material for the note in `allgenus/note/` (key read_D1). Not note text: it carries sources, status labels and project names
that the note itself must not contain (binding rule 5). Every item names its source (file + section) and its status label in the
sources. Pieces: this file (s.0 sources, s.1 definitions); `DIGEST_D1_2.md` (s.2 statements); `DIGEST_D1_3.md` (s.3 proofs);
`DIGEST_D1_4.md` (s.4 conventions, s.5 evidence, s.6 literature, s.7 open issues).

**Source abbreviations** (paths relative to `formalquintic/hodge_vertex_quintic/`; "allgenus/" omitted below `allgenus`):

| abbr. | file |
|---|---|
| MS4, MS5 | `allgenus/MECHANISM_STAGE4.md`, `allgenus/MECHANISM_STAGE5.md` |
| SYN4, SYN4b | `stage4/STAGE4_SYNTHESIS.md`, `stage4/STAGE4_SYNTHESIS_2.md` |
| EVID, ER | `stage4/evid/EVID.md`, `stage4/evid/referee/REFEREE.md` |
| POLY, PR | `stage4/poly/POLY.md`, `stage4/poly/referee/REFEREE.md` |
| GAM, GP, GR | `stage4/gamma/GAMMA.md`, `stage4/gamma/GAMMA_PROOFS.md`, `stage4/gamma/referee/REFEREE.md` |
| PB, PB2 | `stage5/pb/PB_AUDIT.md`, `stage5/pb/PB_AUDIT_2.md` |
| PS, PS2 | `stage5/referees/pb_steps/REFEREE.md`, `.../REFEREE_2.md` (corrections PS-C1..C10 are in PS2 s.6) |
| SYN5, SYN5b | `stage5/STAGE5_SYNTHESIS.md`, `stage5/STAGE5_SYNTHESIS_2.md` |
| RS, RS2, RS3 | `stage5/referees/REFEREES_SYNTHESIS.md`, `_2.md`, `_3.md` |
| UR | `natural_form/theory/UNIFORM_REDUCTION.md` |
| GI | `allgenus/gauge/GAUGE_IDENTITY.md` |
| FSL, FSL2, FSL4 | `stage5/lemmas/FORMAL_SP_LEMMAS.md`, `_2.md`, `_4.md` |
| RL | `stage5/lemmas/referee/REFEREE.md` |
| CP, CP3 | `stage5/chart/CHART_PROOF.md`, `stage5/chart/CHART_PROOF_3.md` |
| RF2, RA, LF | `stage5/chart/referee_formal/REFEREE_2.md`, `stage5/chart/referee_adversary/REFEREE.md`, `stage5/referees/lmi_formal/REFEREE.md` |
| SYN2b | `stage2/STAGE2_SYNTHESIS_2.md` |
| CYC | `stage5/critic/CYC_PROOF.md` (not read by me; see s.0) |

Status labels are those of the sources: PROVED (argument written; "refereed once/twice" = AI referees), VERIFIED (exact
computation, range and field stated), CONJECTURAL, SKETCH. No step has been checked by a human (RS2 s.6.3). Contradiction and
overclaim tags (K1.., O1..) refer to SYN4b s.4-s.5 when prefixed "Stage 4", otherwise to SYN5b s.5-s.6; A1..A10 are RS s.2.

## 0. Sources read, and what I did not read

* **Read in full:** MS4, MS5, SYN4, SYN4b, EVID, ER, POLY, PR, GAM, GP, GR, PB, PB2, PS, PS2, SYN5, SYN5b, RS, RS2, RS3.
* **Read in part** (for the exact definitions and statements the chain cites): UR s.1-s.4; GI s.1-s.3; FSL s.0-s.1.5; FSL2 s.3-s.5
  (to the first lines of s.5.4); FSL4 s.9-s.10; RL s.3.2-s.3.4; CP s.0-s.1.7; CP3 s.3.9-s.3.11; RF2 s.5-s.6 (C1-C5); RA s.3.4;
  LF s.4.2-s.4.3; SYN2b s.6.3; the root `refGW/refs.bib` (by search, for keys only); the head of the conifold note's `LABELS.md`
  (prefix convention).
* **Not read:** CYC and its two referees (cyc_formal, cyc_adversary): known to me only through RS s.2-s.3, RS2 s.5 items 1-16 and
  SYN5b s.9; `stage5/lit/LIT*`; `stage5/general/GENERAL*`; CP s.2-s.3.8, `CHART_PROOF_2.md`, `_4.md`; RA s.0-s.3.3 and RA's part 2;
  RF part 1; LF s.0-s.4.1 and LF part 2; the lmi_adversary reports; FSL s.2, s.6-s.8, s.11-s.12, `FORMAL_SP_LEMMAS_3.md`, `_5.md`;
  RL s.0-s.3.1 and RL part 2; UR s.0, s.5-s.10; GI s.0, s.4-s.9; `GAUGE_IDENTITY_HAT.md`; `NATURAL_FORM.md`; the Stage-2 notes other
  than SYN2b s.6.3. No code and no log was opened or run; every number below is quoted from a note. None of the cited papers was read.
* **Supersession.** MS5 s.2 supersedes the input lists of MS4 s.2 (MS5 header; O10, O11, O18). Within Stage 4, the critic's verdicts
  K1-K11 and corrections O1-O18 (SYN4b) override the worker notes; within Stage 5, SYN5b s.5-s.6 and RS2 s.5 override PB and CP.

## 1. Definitions (final form)

Format: **D1.k name.** definition. *Notation; source.*

**D1.1 Variables and the three points.** `q` the parameter of the formal quintic; `u := 1 - 5^5 q`; `L := u^(-1/5) = (1 - 5^5 q)^(-1/5)`;
`s := u^(1/5) = 1/L`; the hatted variable `w` (also written `zhat`) with `z = s^4 w = u^(4/5) w`; `zeta := e^(2 pi i/5)`.
`D := q d/dq = -(1-u) d/du`; `th_u := u d/du`, `th_w := w d/dw`. Points: large radius `q = 0` (`u = 1`, `L = s = 1`, `z = w`);
conifold `u = 0` (`L = oo`); orbifold `q = oo` (`L = 0`). Colour `i` is the critical point `sigma = zeta^i L`.
*POLY header, s.1; GAM header and s.0 table; PB header; ER s.4; CP s.0.1-s.0.2. UR writes `th := z d/du` (clash, s.4).*

**D1.2 Picard-Fuchs operators.** `M := zD`; `PF5 := M^5 - 1 - q prod_(k=1..5)(5M + kz)`; `PF4 := D^4 - 5q prod_(k=1..4)(5D + k)`;
`PF5 = z^5 D PF4 - 1` (from `q(5D + 5) = 5Dq`). In `u`: `PF5 = u M^5 - 1 - (1-u)(3zM^4 + (17/5)z^2 M^3 + (9/5)z^3 M^2 + (274/625)z^4 M
+ (24/625)z^5)`, coefficients `e_k(1/5, .., 5/5)`. *UR s.1; PB s.2.2; POLY s.1.*

**D1.3 Conifold frame data.** `varpi_0 = C0`: a power-series solution of PF4 in the family `varpi_0 = f_0 + alpha g + beta f_2`
("conifold frame `(alpha, beta)`"; `f_0, f_1, f_2` the Frobenius basis at `u = 0`); `g = f_1 + (7/10) f_2` the vanishing period
(`gamma = 7/10`); `t_c := g/varpi_0 = u + (7/10)u^2 + (41/75)u^3 + (9/20)u^4 + ..` (frame `(0,0)`); `C1 := D t_c = -(1-u) t_c'`. In this
family `varpi_0(0) = 1`, `t_c'(0) = 1`, `C0(0) C1(0) = -1`. `tau(u) := 1/2F1(1, 4/5; 9/5; u)`, `1/tau = Sum_n 4u^n/(5n+4)`;
`v(u) := u tau^(-5/4) = u + (5/9)u^2 + (220/567)u^3 + ..`; `V := tau^(-5/4)` (so `v = uV`; PB also uses `V` for a residue, s.4).
*UR s.1, s.3 Remark (c); GI s.1; PB s.2.4, s.2.8; ER s.1, s.3; RS2 items 1, 27-28 (the family); RS3 s.7.2 (`t_c` to `u^4`).*

**D1.4 The quintic cyclic system.** `zD Y_a = c_a Y_(a+1)` (`a` in Z/5), `c_a = (C1, L^5/(C0^2 C1^2), C1, C0, C0)`; equivalently
`th Y = A Y` (`th = z d/du`), `A = Sum_a b_a E_(a,a+1)`, `b_a = -c_a/(1-u)`, i.e.
`b = (t_c', -1/(u varpi_0^2 (1-u)^3 t_c'^2), t_c', -varpi_0/(1-u), -varpi_0/(1-u))`: units except `b_1`, with `u b_1` a unit ("one
simple pole"); `Prod b = -1/(u(1-u)^5)`. Its closing relation `4 -> 0` is PF5 by (CYC) (other digest). *UR s.1; GI s.1; PB s.2.3, s.2.8.*

**D1.5 The model and its WKB rows.** `lambda := -(x^4/4 + v/x)`, `I_m := Int x^m e^(lambda/z) dx`; rows `Y = (I_0, -I_(-1), -I_3, I_2, I_1)`
solve `z d/dv Y_b = bmod_b Y_(b+1)`, `bmod = (1, -1/v, 1, -1, -1)`. Exponents `c_b = (0, -1, 3, 2, 1)` (`b = 0..4`; `c_b = -b mod 5`),
signs `eps_b = (+, -, -, +, +)`, so `Y_b = eps_b I_(c_b)`. Prop. S: `Y_b = e^(U/z) sigma^(c_b - 1) F_b(t)`, `t = z/sigma^4`,
`U = -(5/4) zeta^i sigma^4`, `sigma = v^(1/5)`; T-system `F_1 = T_(-1) F_0`, `F_2 = -T_(-2) F_1`, `F_3 = T_2 F_2`, `F_4 = -T_1 F_3`,
`F_0 = -T_0 F_4`, `T_k = -zeta^i + (t/5)(k - 4t d/dt)`; unique solution with `F_0(0) = 1`; `F_b(0) = eps_b`;
`F_0 = 1 - (3/20)t + (9/800)t^2 + (11571/80000)t^3 - ..`; `F_(0,n) := [t^n] F_0`. At `v(u)`: `Q^m_b = sigma^(c_b - 1) F_b(z/sigma^4)` with
`sigma = s tau^(-1/4)`, `z/sigma^4 = tau w`. Model matrix `R^(t) = Psi0 M_0(t)`, `R^_b^(i)(t) = zeta^(i(b+1)) F_b(zeta^-i t)`.
*UR s.1-s.2, s.4; PS s.1 (ii)-(vi); PB s.2.6; EVID s.2.2.*

**D1.6 Vanishing solutions.** `I^C_m(v, z) := Sum_(i>=0) (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!` (terms with `m + 4i + 1 >= 0`)
`= Res_(x=0) x^m e^(-(x^4/4 + v/x)/z) dx`; `Y^C := (I_0, -I_(-1), -I_3, I_2, I_1)^C`, i.e. `Y^C_b = eps_b I^C_(c_b)`; `Y^C(0, z) = -e_1`;
`Y^C(v(u), z)` in `Q[z^-1][[u]]^5`. Quintic: `y^(0) := -e_1`, `y^(n+1)_a := Int_0^u b_a y^(n)_(a+1) du`, `Y^C_Q := Sum_n z^-n y^(n)` in
`Q[z^-1][[u]]^5`; `y^(n)_a = 0` unless `n = 1 - a mod 5`, and then `ord_u y^(n)_a = (4n + c_a + 1)/5`; `[z^-1] Y^C_Q = -t_c e_0`.
*GI s.2 (Lemma V, Lemma Q).*

**D1.7 Regular gauge, normalisations, (R1).** With `At := v'(u) Amod(v(u))` (the pulled-back model): `P = Sum_k z^k P_k(u)`,
`P_k` in `gl_5(Q[[u]])`, `P_0` diagonal with unit entries, `z P' = A P - P At`; graded: `P_k = Sum_b p_(k,b) E_(b+k,b)`; exists by
Theorem U, unique up to `P -> P c(z^5)`. Canonical `P_can`: `P_11(0, z) = 1`. Stored: `p_(5q,0)(0) = 0` (`q >= 1`); then
`P^s = P_can (1 + m5 z^5 + ..)`, `m5 = -43397/807975000`. (R1): column 1 of `P(0, z)` is `P_11(0, z) e_1`, row 2 is `P_22(0, z) e_2^T`, and
`P_22(0, z) = P_11(0, z)`. Hatted: `Ghat_ab = tau^((1-c_b)/4) Sum_k w^k u^(mu(b,k)) p_(k,b)(u) delta_(a,b+k)`,
`mu(b, k) = (4k + c_b - c_(b+k))/5` (an integer `>= 0`). *UR s.3-s.4; GI s.1, s.3 (with Corollary (a)); EVID s.1, s.5.*

**D1.8 The normalisation constant `D`; (F_g), (H), (N).** `Q` := the actual (QRR-normalised) WKB fundamental solution of the quintic
system. R0: `Q D = P_can Q^m(v(u))` with `D = diag(D_i)` constant in `u`, in `Q[[z]]`, defined after stripping the SAME exponential
on both sides (UR's `e^((u_i - u_i^c)/z)` differs from `e^(u_i/z)` by a constant that is not a power series in `z`), and
`D_i(z) = D_0(zeta^-i z)`. **(F_g)**: `D_0 = 1 mod z^(3g-2)` (what Conjecture N(g) needs). **(H)**: `D_0` in `Q[[z^5]]`. **(N)**: `D_0 = 1`.
*UR s.4; SYN4 header; MS4 s.1; PS s.4.1; RS2 item 24. UR s.4's "(H) purity" is the same statement read on the rows (clash, s.4).*

**D1.9 The colour-0 WKB function.** `Phi(w, u) = Sum_n w^n ph_n(u)`, `ph_0 = 1` (the hatted colour-0 entry of row 0, times `C0`);
flat amplitude `A_n := L^(4n+1) ph_n(L^-5) = u^(-(4n+1)/5) ph_n` (`psi_0 = e^(u_0/z) Sum_n A_n z^n`, `D u_0 = L`); `R_n(L) := L^(4n) ph_n(L^-5)
= A_n/L`. Colour `i`: `Phi^(i)(w, u) = Phi^(0)(w/zeta^i, u)`. HPF: the hatted form of PF5 (s.2, `lem:pol-HPF`).
*GAM s.1; POLY header, s.1; PB s.3; PS s.4.2 (a).*

**D1.10 Normalisations of `Phi`; free coefficients.** `f_m := [u^(4m)] ph_(5m)` (EVID writes `f_(5m)`, GAM `x_m`). `Phi^(0)`: all
`f_m = 0`. `Phi^QRR`: `A^QRR_n(q = 0) = h_n` for every `n` (PB's H9; equivalently `Phi^QRR(z, 1) = h(z)`); EVID fixes `f_(5m)` by
`ph_(5m)(1) = h_(5m)` (EVID s.2.1 (c)). `Phi^can := C0 x` (colour-0 entry of row 0 of the hatted `P_can Q^m(v(u))`) (PB s.2.1);
equivalently the polynomial solution with `K_q = 0` for all `q` (PS s.4.4). `Phi^SP`: D1.13. `F_inf(z) := Sum_m f_m z^(5m)`, the value of
the colour-0 WKB at `L = 0`. *POLY s.2-s.3; EVID s.2.1, s.2.3; PB s.1, s.2.1; PS s.4.4; ER s.4.*

**D1.11 QRR series; `Ginf`.** `b_k := B_(2k)/(2k(2k-1))`; `h_w(z) := exp(-Sum_k b_k z^(2k-1) Sum_(a in w) a^(1-2k))`;
`w_i := {zeta^i - zeta^j}_(j != i) U {-5 zeta^i}`; `h := h_(w_0)`, `h_n := [z^n] h`; equivalently
`h(z) = exp(Sum_k b_k z^(2k-1) [5^(1-2k) + Sum_(j>=1) (zeta^j - 1)^(1-2k)])`; `h_1 = -3/20`, `h_5 = -257071/8064000000`;
`h_(w_i)(z) = h_(w_0)(zeta^-i z)`. `Ginf(z) := exp(Sum_(l>=1) (-1)^(l+1) B_(5l+1)(1/5) z^(5l)/(l(5l+1)))`, the normalised Stirling
factor of `prod_j Gamma(zeta^j/z + 1/5)`. *MS4 s.1, s.4.3; MS5 s.1-s.2; ER s.6; POLY s.1; GAM s.0 item 6; PS s.4.5; LF s.4.2.*

**D1.12 The flat full-Stirling saddle series `c(z; L)`.** For `L` in `Lam := C \ {0, 1, zeta, .., zeta^4}`, `q(L) := (1 - L^-5)/5^5`. On a
disc around `L`: `lambda := -5 log(5L) + Sum_j log(zeta^j - L)` (`e^lambda = -q(L)`); `S_0(sigma) := lambda sigma + 5 sigma log(5 sigma) -
5 sigma + Sum_j [(zeta^j - sigma) log(zeta^j - sigma) - (zeta^j - sigma)]`; `S_1 := (1/2) log(5 sigma) - (1/2) Sum_j log(zeta^j - sigma)`;
`S_2k := b_k [(5 sigma)^(1-2k) + Sum_j (zeta^j - sigma)^(1-2k)]`; `A_St := exp(S_1 + Sum_k z^(2k-1) S_2k)`.
`c(z; L) := Nhat[S_0, A_St; L](z)`, the normalised formal SP at `sigma = L` (measure `d sigma`), `c_0 = 1`. It is the normalised saddle
series of `Gamma(1 + 5s) prod_(j=0..4) Gamma(zeta^j/z - s) (-q)^s d sigma` (`sigma = zs`) with all six Gamma factors replaced by `Gamma_St`
(the `(log z)/z` terms cancel since `Sum_j zeta^j = 0`). `S_0'' = 5/(sigma(1 - sigma^5))`; saddles `u sigma^5 = 1`.
*CP s.0.1, s.1.6-s.1.7; GAM s.2; GR s.2.2.*

**D1.13 The 1-form, `Phi^SP`, and the residue `V`.** In `x = 1/sigma`: phase `W = Lam/x + G(x)`, `Lam := -log(1-u)`,
`G(x) := Sum_j (1/x - zeta^j) log(1 - zeta^j x) = Sum_(m>=1) x^(5m-1)/(m(5m-1))`; amplitude `at(x, z)`, `log at = -(1/2) log(1 - x^5) +
Sum_k b_k (xz)^(2k-1) [5^(1-2k) - 5 Sum_(P>=0) C(2k-2+5P, 5P) x^(5P)]`, `at(0, z) = 1`. Hatted (`x = sX`, `z = s^4 w`): `W/z = What/w`,
`What(X; u) = Sum_(m>=1) (u^(m-1)/m)[1/X + X^(5m-1)/(5m-1)]`, `B(X) := at(sX, s^4 w)` (polynomial in `X`, `w` at each `u`-order;
`B(1)|_(w=0) = (1-u)^(-1/2)`); `What'(X) = -(uX^2)^(-1) log((1-uX^5)/(1-u))`, critical point `X = 1` for every `u`,
`What''(1) = 5/(1-u)`, `What(1) = (5/4)/tau`. **`Phi^SP := SP_(X=1)[e^(-What/w) B dX] / (sqrt(2 pi w/5) e^(-What(1)/w))`**, `ph^SP_0 = 1`.
**`V := Res_(x=0) e^(-W/z) at dx = Sum_(k>=1) (-Lam/z)^k/k! [x^(k-1)](e^(-G/z) at)`**, in `Q[z^-1][[u]]`. *PB s.2.1; POLY s.5; PS s.2.2.*

**D1.14 CFU map, reduction coefficients, `P^SP`.** CFU: `Y = X(1 + Sum_(k>=1) u^k Y_k(X^5))`, `Y_k` in `Q[T]`, `T = X^5`, with
`What(X; u) = f(Y) := V/Y + Y^4/4`, `V = tau^(-5/4)`; `Y_c := Y(1; u) = V^(1/5) = tau^(-1/4)`. Reduction: `beta(Y) dY := B(X(Y)) X'(Y) dY`,
`beta` in `Q[Y, w][[u]]`; `beta e^(-f/w) dY = Sum_(m=-1..3) pi_m Y^m e^(-f/w) dY + d(Gamma e^(-f/w))`, `pi_m` in `Q[w][[u]]`, `Gamma` in
`Y Q[Y, w][[u]]`. `Q^SP`: rows from `Phi^SP` by the cyclic system, colours by covariance; `P^SP := Q^SP (Q^m)^-1`; row-0 entries
`G_b := s^(c_b) P_0b`. *PB s.2.1, s.2.4-s.2.6; POLY s.5-s.6.*

**D1.15 The formal functionals** (defined in the formal-SP appendix, another digest). `SP^norm_c[F; a] := E_(zC)[exp(F_3/z) a(c + t)]`,
`C = -H^-1` (FSL Def. 1.2); full SP with the symbols `[e^(F(c)/z)] (2 pi z)^(n/2) [Delta^(-1/2)]` (FSL Def. 1.3). PB's one-variable
definition through the Morse coordinate, `SP[g dt] := sqrt(2 pi w/a) e^(-f_0/w) Sum_m (2m-1)!! (w/a)^m [y^(2m)](g(t(y)) t'(y))`, is the same
functional (PB s.2.0, FL2; FSL2 s.5 Remark). Phase conventions: `e^(Phi/z)` (GAM) is `F = Phi`; `e^(-W/w)` (POLY) is `F = -W`, `z = w`
(FSL s.1.4). Residue: ring `L := {Sum_(j in Z) c_j X^j : c_j in K((v)), ord_v(c_j) >= eps|j| - M}`, `v = 1/w`; `Res := [X^-1]`,
coefficientwise in `u` (FSL4 s.10.1). *FSL s.1.4; FSL4 s.10.1; PB s.2.0, s.2.7; PS s.2.1 (PS-C2).*

**D1.16 Objects of the finite route.** For any polynomial-type colour-0 `Phi` with hatted rows `Qhat` (all colours by covariance):
`Ghat := Qhat R^(tau w)^-1`, `K_q := [u^(4q) w^(5q)] Ghat_11`. `S_q(u) := Sum_(n<=5q) (A_(-1) Phi)_n (-tau)^(5q-n) F_(2,5q-n)`,
`(A_(-1) Phi)_n = ph_n + (1-u)((4n-3)/5 - th_u) ph_(n-1)`. First row: `Phi/C0 = Sum_k w^k g_k(u) F_(-k mod 5)(tau w)`,
`g_k = u^(mu(-k,k)) tau^((1-c_b)/4) p_(k,b)(u)`, `b = -k mod 5`, `mu(-k, k) = (4k + c_(-k))/5 = 4q + (0,1,2,3,3)_r` for `k = 5q + r`.
*EVID s.2.2-s.2.3; ER s.3; POLY s.3 (the `mu` table corrected there from a first draft `4q + r`).*
