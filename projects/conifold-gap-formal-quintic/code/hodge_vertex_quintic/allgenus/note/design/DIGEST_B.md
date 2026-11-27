# DIGEST_B (key read_B): the two cyclic systems, Theorem U, the gauge identity, symplecticity, (CYC)

Design-phase digest for the note in `allgenus/note/`. Author-facing only; no LaTeX. Subject: the quintic's cyclic system (rows,
coefficients `c_a`, PF4, PF5), the model system of `x^4 + v/x` (rows, `v(u)`, UR Prop. S), Theorem U with its proof, the canonical
gauge `P_can` and `P_11(0,z) = 1`, `Q D = P_can Q^m(v(u))` and `D_i(z) = D_0(zeta^-i z)`, GI Theorem G with Lemmas 1, V, Q,
symplecticity, the residue relations, and (CYC) with its proof.
Pieces: `DIGEST_B.md` (s.0, s.1), `DIGEST_B_2.md` (s.2 S1-S19), `DIGEST_B_3.md` (s.2 S20-S38), `DIGEST_B_4.md` (s.3),
`DIGEST_B_5.md` (s.4, s.5), `DIGEST_B_6.md` (s.6, s.7).

**Abbreviations** (paths relative to HVQ = `formalquintic/hodge_vertex_quintic/`).
UR = `natural_form/theory/UNIFORM_REDUCTION.md`; NFS = `natural_form/SYNTHESIS.md` (critic of the natural-form round; the only written
audit of UR); RN = `natural_form/referee/REFEREE_NAT.md`; VT = `natural_form/verify_theory/` (logs only, no report); LG =
`natural_form/lg/LG_MODEL.md`; GI = `allgenus/gauge/GAUGE_IDENTITY.md`; GH = `allgenus/gauge/GAUGE_IDENTITY_HAT.md`; GR =
`allgenus/gauge/referee/REFEREE.md`; MEC = `allgenus/MECHANISM.md`; S01 = `allgenus/SYNTHESIS_STAGE01.md`; MS2, MS4, MS5 =
`allgenus/MECHANISM_STAGE2.md`, `_STAGE4.md`, `_STAGE5.md`; ST2, ST2b = `allgenus/stage2/STAGE2_SYNTHESIS.md`, `_2.md`; ST4b =
`allgenus/stage4/STAGE4_SYNTHESIS_2.md`; ST5, ST5b = `allgenus/stage5/STAGE5_SYNTHESIS.md`, `_2.md`; CYC =
`allgenus/stage5/critic/CYC_PROOF.md`; CF, CF2, CF3 = `allgenus/stage5/referees/cyc_formal/REFEREE.md`, `_2.md`, `_3.md`; CA, CA2 =
`allgenus/stage5/referees/cyc_adversary/REFEREE.md`, `_2.md`; RS1, RS2 = `allgenus/stage5/referees/REFEREES_SYNTHESIS.md`, `_2.md`;
PB, PB2 = `allgenus/stage5/pb/PB_AUDIT.md`, `_2.md`; PS = `allgenus/stage5/referees/pb_steps/REFEREE.md`; EV = `allgenus/stage4/evid/EVID.md`;
TH = `allgenus/theory/REDUCTION_DRAFT.md`; FQ2 = `FORMAL_QUINTIC_genus2.md`; CS = `CONIFOLD_SYMBOL_genus2.md`.
**Status labels** are the sources': PROVED, PROVED-mod(X) (NFS), VERIFIED (exact finite computation, range stated), QUOTED, SKETCH,
CONJECTURAL, REFUTED. "AI-audited/AI-refereed" is recorded here only so the writers know the review state; rule (5) keeps it out of
the note. Final status = after MS5 s.2 and RS2 s.6 (MS5 supersedes MS4 s.2; MS4 updates MS2 s.1; MS2 supersedes MEC s.1, s.3, s.4).

## s.0 Sources read, and what was not read

**Read in full.** UR s.0-s.10 and "Sources (web)". NFS s.0-s.6 (item 2 records the audit of Theorem U: "audited line by line by a
second AI and re-implemented along the proof route"). RN s.0-s.7 (written before Theorem U existed, NFS s.2 item 2). GI s.0-s.9.
GH s.0-s.6. GR s.0-s.5. MEC s.0-s.6. S01 s.0-s.7. MS2 s.1-s.6. MS4 s.1-s.5. MS5 s.1-s.5. CYC s.0-s.9. CF s.0-s.5, CF2 s.6-s.9,
CF3 s.10-s.12. CA s.0-s.3, CA2 s.4-s.8. RS1 s.0-s.4. RS2 s.5-s.6. ST4b s.4-s.7. ST5b s.5-s.10.
**Read in part.** ST2 s.2 (the chain and the input table); ST2b s.4-s.6 (contradictions K1-K11, overclaims O1-O27, write-up plan
s.6.3); ST5 s.1-s.2, s.4.1 (input table, chain S1-S9, the critic's (CYC) checks); PB s.0-s.3 (s.2.3 (CYC) outline, s.2.8 hypotheses
of U, s.3 the glue); PB2 s.4-s.5 (code table, findings P1-P5); PS s.3-s.4 (s.2.8-s.2.9 and the glue); EV s.2.1-s.2.4 (Lemma N);
TH s.0-s.1.5 (hypotheses (Sp), (A) on a dressing; the canonical gauge as example (i)); LG s.0-s.3 (model conventions); FQ2 s.1-s.3;
CS s.1 (first half). Logs: `natural_form/theory/ur_model_check.log`, `ur_chain_check.log`, `ur_gauge_run.log`; VT
`v1_model.log`, `v2_run.log`, `v3_fourier.log`, `v4_symplectic.log` and `v2_run.py` lines 1-10, 103-124 (its range K = 8);
`natural_form/critic/critic_gauge_vs_table.log`. Templates (format only): `resolvedconifold/hodge_vertex_conifold/note/OUTLINE.md`
s.0, `NOTATION.md` (head), `LABELS.md` (head), `DRAFT_NOTES.md` (head).
**Not read** (outside this key, or not needed for it): `OU_EXTENSION_genus3.md` and `NATURAL_FORM.md` (hatted conventions: the matrix
`Psi(u)`, the kinematic plane, the outer graph sum; items depending on them are marked [conv: OU3 s.1]); `stage4/poly/POLY.md` (its
Lemma 2 is quoted here from PB s.1 H7 and PS s.3); `theory/REDUCTION_DRAFT_2.md`, `_3.md` (transport Lemmas E, T, F, Theorem R);
`univ/`, `model/`, `tr/`, `xyswap/`, `stage2/*` worker notes, `stage4/gamma/`, `stage5/chart/`, `lemmas/`, `general/`, `lit/`;
RS3 (only grepped: its s.8 human-check list); VT has no written report (NFS s.0 summarises it); the code (`fq_*.py`, `ur_*.py`,
`gi_*.py`, `cyc_*.py`) apart from the lines named above; no paper was opened (s.6 records what the sources say about each).

## s.1 Definitions (final form; names are proposals; clashes in s.4.3)

**D1 Variables.** `q` near 0 (large radius), `u = 1 - 5^5 q` (conifold `u = 0`, large radius `u = 1`), `s = u^(1/5) = 1/L`,
`L^5 = 1/u = 1/(1 - 5^5 q)`, `D = q d/dq = -(1-u) d/du`, `z` the descendant variable (a constant for `D`), `th = z d/du`, so
`zD = -(1-u) th`; hatted `zhat = z/s^4 = L^4 z` (PB, POLY write `w` for `zhat`); `zeta` a primitive fifth root of unity; colour `i` in Z/5.
Source: UR header, s.1; CYC s.0; CS s.1; PB header.

**D2 PF4 and PF5.** `PF4 := D^4 - 5q prod_(k=1..4)(5D + k) = Sum_(k=0..4) A_k D^k`, `A_4 = u`, `A_3 = -2(1-u)`, `A_2 = -(7/5)(1-u)`,
`A_1 = -(2/5)(1-u)`, `A_0 = -(24/625)(1-u)`. `PF5 := M^5 - 1 - q prod_(k=1..5)(5M + kz)`, `M = zD`; u-form
`PF5 = u M^5 - 1 - (1-u)(3zM^4 + (17/5)z^2M^3 + (9/5)z^3M^2 + (274/625)z^4M + (24/625)z^5)` (coefficients `5^(5-j) e_j(1..5)/3125`);
`PF5 = z^5 D PF4 - 1` (S1). PF5 annihilates the I-function of the formal quintic at the zeta point (FQ2 s.1; CA s.2.7 re-derives
the ratio of the d-th and (d-1)-th terms). Source: CYC s.0 and Cor.; CF s.1 (A1); CA s.2.4, s.2.7; UR s.1 (u-form with `M = -(1-u)th`).

**D3 The ring and the solution space at u = 0.** `K := Q((u))` with `D`; `R := K[l]`, `l` an indeterminate, `D l := 1 - 1/u` (`l = log u`).
Power-series solutions `f_0, f_1, f_2`, `f_k = u^k + O(u^3)`: `f_0 = 1 + (2/625)u^3 + ..`, `f_1 = u - (19/60)u^3 + ..`,
`f_2 = u^2 + (37/30)u^3 + ..`; log solution `y_3 = g l + h`, `h` in `Q[[u]]` with `h_0 = h_1 = h_2 = 0`
(`h = -(23/360)u^3 - (6397/60000)u^4 + ..`), `g := f_1 + gamma f_2`, `gamma = 7/10`: the **vanishing period**
`g = u + (7/10)u^2 + (41/75)u^3 + (1133/2500)u^4 + ..`. `S := span_Q(f_0, f_1, f_2, y_3)` (the solutions used; CF2 s.6: no claim that
S is the whole solution space in R is needed). Source: CYC s.0-s.1; CF s.2.2-s.2.3; CA s.1; FQ2 s.1.

**D4 Conifold frames.** Programme family (fq_conifold): `varpi_0 = f_0 + alpha g + beta f_2`, `(alpha, beta)` in Q^2, `t_c = g/varpi_0`;
then `varpi_0(0) = 1`, `t_c = u + gamma_t u^2 + ..`, `t_c'(0) = 1`, `D t_c(0) = -1`, `w_1 := [u^1] varpi_0 = alpha`, `gamma_t + w_1 =
gamma_g = 7/10`. Default frame `(0,0)`: `varpi_0 = 1 + O(u^3)`, `t_c = u + (7/10)u^2 + (41/75)u^3 + (9/20)u^4 + ..`.
Wider CYC frame (RS1 A1): `y_0` any non-zero power-series solution of PF4 not proportional to `g`, `t = g/y_0`.
Label warning (RS1 A2 = CF-E3 = CA-K1): the critic's checks use `y_0 = f_0 + a f_1 + b f_2`; critic `(1/3,-2/7)` = programme
`(1/3,-109/210)`, critic `(7/11,-3/5)` = programme `(7/11,-23/22)`; programme `(1/3,-2/7)` = `f_0 + f_1/3 - (11/210) f_2`.
Source: CYC s.0; CA s.1, CA2 s.4; CF2 s.9; UR s.1, s.5; GI s.1; RS2 s.5 items 3, 5, 28.

**D5 The quintic cyclic system.** `c = (c_0, .., c_4) := (Dt, L^5/(varpi_0^2 (Dt)^2), Dt, varpi_0, varpi_0)` with `t = t_c`; `prod c_a = L^5`.
System `zD Y_a = c_a Y_(a+1)`, `a` in Z/5; equivalently `th Y = A Y`, `A = Sum_a b_a E_(a,a+1)`, `b_a = -c_a/(1-u)`:
`b = (t_c', -1/(u varpi_0^2 (1-u)^3 t_c'^2), t_c', -varpi_0/(1-u), -varpi_0/(1-u))`, `prod b_a = -1/(u (1-u)^5)`;
`A = u^-1 r E_12 + A_0 + u A_1 + ..`, `r = -1`, `A_j` constant rational matrices. Self-dual: `b_0 = b_2`, `b_3 = b_4`, i.e. `eta A = A^T eta`.
*Rows.* `fq_rows`: `X_0 = psi_i/C_0`, `X_(a+1) = (1/(C_1 kappa_a))(zD + zeta^i L) X_a` (`a = 0..3`), `kappa = (1, L^5/(C_0^2 C_1^3), 1, C_0/C_1,
C_0/C_1)`, where `e^(u_i/z) z psi_i` is the colour-i WKB solution of PF5 (`D u_i = zeta^i L`); `Y_a = e^(u_i/z) X_a`; in the conifold
frame `C_0 -> varpi_0`, `C_1 -> D t_c = -(1-u) t_c'`; hence `1/kinv[a] = C_1 kappa_a = c_a` (`a = 0..3`). The closing coefficient
`c_4 = C_1 kappa_4 = C_0` is used by no code; it is forced by the symbol: `prod c_a = L^5` (CA s.2.1). The normalised rows
`Q_a^(i) = zeta^i X_a` satisfy the same system (UR s.1). Source: UR s.1; CYC s.0 Cor.; CA s.2.1-s.2.2; CF2 s.8.2; GI s.1.

**D6 Canonical coordinates.** Eigenvalues of `A`: `du_i/du = -zeta^i u^(-1/5)/(1-u)`; `u_i - u_i^c = -zeta^i (5/4) u^(4/5) 2F1(1,4/5;9/5;u)
= -zeta^i (5/4) u^(4/5)/tau(u)`. Source: UR s.1; RN s.3(a) (`Int_0^u t^(-1/5)(1-t)^(-1) dt = (5/4) u^(4/5) 2F1(1,4/5;9/5;u)`).

**D7 tau and v.** `tau(u) := 1/2F1(1,4/5;9/5;u) = 1/Sum_(n>=0) 4u^n/(5n+4) = 1 - (4/9)u - (50/567)u^2 - (4300/96957)u^3 - ..`;
`v(u) := u tau^(-5/4) = u + (5/9)u^2 + (220/567)u^3 + (29030/96957)u^4 + ..`: the unique power series with `v = u + O(u^2)` and
`v'^5/v = 1/(u(1-u)^5)` (equivalently `v^(-1/5) dv = u^(-1/5) du/(1-u)`, `v^(4/5) = u^(4/5)/tau`); without the normalisation `v'(0) = 1`
only up to `v -> c v`, `c^4 = 1` (NFS s.4 "only up to a fourth root of unity"; PS s.3 "picks the fifth root"). `sigma := v^(1/5) =
s tau^(-1/4)`, so `z/sigma^4 = tau zhat`. Source: UR s.0 item 4, s.3 Rmk (c), s.4; PS s.3; `ur_gauge_run.log` (series).

**D8 The model and its system.** LG model `lambda(x) = -(x^4/4 + v/x)`, primitive form `dx`; formal integrals `I_m(v,z) = Int x^m
e^(lambda/z) dx` with (Mod1) `z d/dv I_m = -I_(m-1)` and (Mod2) `m z I_(m-1) - I_(m+3) + v I_(m-2) = 0`. Rows
`(Y_0, .., Y_4) := (I_0, -I_(-1), -v I_(-2) = -I_3, I_2, I_1)`; model system `z d/dv Y_b = bmod_b Y_(b+1)`, `bmod = (1, -1/v, 1, -1, -1)`;
scalar equation `(v th_v^5 + 3 z th_v^4 + 1) Y_0 = 0`, `th_v = z d/dv`. Dictionary of classes: `(1, H, H^2, H^3, H^4) <-> (dx, -dx/x,
-x^3 dx = -v x^(-2) dx, x^2 dx, x dx)`, i.e. `H^a <-> N_a x^(c_a) dx`, `c = (0,-1,3,2,1)`, `N = (1,-1,-1,1,1)`, stated with exponents in the
window -1..3 (not `x^a dx`; "`x^(-a) dx` mod 5" is not an identity of classes). Homogeneity `x -> kx, v -> k^5 v, z -> k^4 z`; critical
points `x^5 = v`. Pulled back by `v(u)`: `At := v'(u) Amod(v(u))`, `bt = (v', -v'/v, v', -v', -v')`, `rt = -1`, `prod bt = -v'^5/v`.
LG writes `lambda = x^4 + u_LG/x`, `u_LG = 4u`, `z_LG = -4z`: the same statement (NFS s.2 item 7). Source: UR s.0 item 2, s.1 "Model";
GI s.1; NFS s.1 item 5, s.2 item 1; LG s.0 item 2, s.2-s.3; PB s.2.8.

**D9 Hatted exponents and weights.** `c = (0,-1,3,2,1)`: `c_(a+1) = c_a - 1 + 5 delta_(a,1)`, `c_a + c_(3-a) = 2`, `c_a = -a + 5 e_a`,
`e = (0,0,1,1,1)`. Hatted rows `Qhat_a = s^(1-c_a) Q_a`. `mu(b,k) := (4k + c_b - c_(b+k))/5` in Z_(>=0), `mu(b, 5q + r) = 4q + mu(b, r)`,
`mu(b,k) >= (4k - 4)/5`. Source: UR header, s.4, s.5 table; GI Lemma Q; GH s.1.

**D10 Cyclic system with one simple pole.** `O = C[[u]]` (UR; also convergent series) or `Q[[u]]` (GI), `K = O[1/u]`: `th Y = A Y`,
`A = Sum_a b_a E_(a,a+1)`, `b_a` in `O^x` for `a != 1`, `u b_1` in `O^x`. Self-dual: `b_0 = b_2`, `b_3 = b_4`. Source: UR s.3; GI s.1;
TH s.1.3 ("any self-dual cyclic one-pole system").

**D11 Regular gauge.** For two such systems `A`, `At`: `P(u,z) = Sum_(k>=0) z^k P_k(u)`, `P_k` in `gl_5(O)`, `P_0` invertible, with
`z P' = A P - P At` (`' = d/du`), so that `Y = P Yt` maps solutions of `th Yt = At Yt` to solutions of `th Y = A Y`. Graded:
`P_k = Sum_b p_(k,b) E_(b+k,b)`. P lies in `gl_5(Q[[z]][[u]])` inside `gl_5(L[[u]])`, `L := Q((z))`. Source: UR s.3; GI s.1.

**D12 Canonical gauge.** `P_can` := the graded regular gauge from `At` (pulled-back model) to `A` (quintic) with `P_11(0,z) = 1`;
equivalently the unique regular gauge with `P Y^C(v(u)) = Y^C_Q` (S28(a)). Code normalisation "stored": `p_(5q,0)(0) = 0`, `q >= 1`;
`P_stored = P_can (1 + m5 z^5 + ..)`, `m5 = -43397/807975000` (unified sign, S01 C11, MEC s.6). Source: UR s.4; GI s.1, s.3 Cor.(a).

**D13 Vanishing solutions.** `I^C_m(v,z) := Res_(x=0) x^m e^(-(x^4/4 + v/x)/z) dx = Sum_(i>=0) (-1/(4z))^i/i! (-v/z)^(m+4i+1)/(m+4i+1)!`
(terms with `m + 4i + 1 >= 0`); `Y^C := (I^C_0, -I^C_(-1), -I^C_3, I^C_2, I^C_1)`, `Y^C(0,z) = -e_1`. Quintic: `y^(0) := -e_1`,
`y^(n+1)_a := Int_0^u b_a y^(n)_(a+1) du`, `Y^C_Q := Sum_n z^(-n) y^(n)`. Both are characterised as the regular solution with initial
value `-e_1`, not by "only negative powers of z" (every regular solution with z-independent initial value is a series in `1/z`; GR s.2).
Source: GI s.2; GR s.2.

**D14 The constant diagonal series D.** `Q` := the matrix of the quintic's WKB columns (rows `a`, colours `i`), normalised by the QRR
constants at `q = 0`; `Q^m(v)` := the model's homogeneous WKB columns (S18); `Q D(z) = P_can Q^m(v(u))` defines `D = diag(D_i(z))`,
constant in u, `D_i` in `Q[[z]]`, with the SAME exponential stripped on both sides (`e^(u_i^c/z)` is not a power series in z and is
not absorbed into D). Source: UR s.4; PS s.4.1; ST2 s.2 (F_g); MS4 s.1.

**D15 Hatted gauge.** `Ghat(u, zhat) := Lambda P(u, s^4 zhat) Lambda'^(-1)`, `Lambda = diag(s^(4-c_a))`, `Lambda' = diag(s^3 sigma^(1-c_b))`;
`[zhat^k] Ghat_(b+k,b) = tau^((1-c_b)/4) u^(mu(b,k)) p_(k,b)(u)`. `G := Psi(u)^(-1) Ghat Psi(0)`, `R' = Psi(u) G`, rows
`Q'(u,x) = Ghat(u,x) Psi0`, `Psi0[b][i] = eps_b zeta^(i(b+1))`, `eps = (1,-1,-1,1,1)`, `Psi0^T eta^(-1) Psi0 = 1` [conv: OU3 s.1 for
`Psi(u)`]. Source: UR s.4, s.6.1; GH s.1.

**D16 Pairing.** `eta_ab = 5 delta_(a+b = 3 mod 5)`, `eta^(-1) = eta/25`; `P(-z)^T eta P(z)` and `P(z) eta^(-1) P(-z)^T` give equivalent
symplecticity conditions. Source: GI s.0-s.1; FQ2 s.1; GR s.2.

**D17 Objects of (CYC).** For `L = Sum_(k<=4) A_k D^k`: formal adjoint `L* := Sum_k (-D)^k o A_k`; concomitant
`P(f,g) := Sum_(k=1..4) Sum_(j=0..k-1) (-1)^j D^j(A_k g) D^(k-1-j) f`; `B(f,g) := P(f,g)` on S; jet matrix `P_ab` with
`B(f,g) = Sum D^a f P_ab D^b g`, for PF4 `P = [[0,(2/5)(1-u),1-u,-u],[-(2/5)(1-u),0,u,0],[-(1-u),-u,0,0],[u,0,0,0]]`;
`Bmat := (B(y_i,y_j))`, `Wr := (D^a y_i)`; the automorphism of R `l -> l + 1` (CYC writes `tau`; CF3 E6 suggests `sigma`);
`th := (Dt)^(-1) D`; `T_k := Sum (Bmat^(-1))_ij y_i D^k y_j`; `M_cyc(y_0, y_1) := (1/y_0) D (1/Dt) D (u y_0^2 (Dt)^2) D (1/Dt) D (1/y_0)`,
`t = y_1/y_0` (CYC writes `M`). Source: CYC s.0-s.6; CF s.4.3; CA s.1; CF3 E6.

**D18 Model WKB rows and the model R-matrix.** Colour i: `Y_b = e^(U/z) sigma^(c_b - 1) F_b(t)`, `t = z/sigma^4`, `U = -(5/4) zeta^i sigma^4`;
`T_k := -zeta^i + (t/5)(k - 4t d/dt)`; model rows at `v(u)`: `Q^m_b = sigma^(c_b-1) R^_b(z/sigma^4)`; `M_0` the model's R-matrix in
Teleman's colour-colour form (LG s.0 item 4, s.2 (D)). Source: UR s.2, s.4; LG s.0, s.2.

**D19 Conifold symbol.** `val_s Q^(i)_(a,n) = -4n + c_a - 1`; `R^_(a,n)^(i) := [s^(-4n+c_a-1)] Q^(i)_(a,n)`, `R^_a^(i)(zhat) = Sum_n
R^_(a,n)^(i) zhat^n`; covariance `R^_a^(i)(zhat) = zeta^(i(a+1)) R^_a^(0)(zeta^(-i) zhat)`. Source: CS s.1.

**D20 Euler-dilaton vector.** `e(z) = z 1 + u H` (flat frame); for the homogeneous model `tau_1(1) + v tau_0(H)` (= TH's `frame(1,0)`).
Convention of record: `H = d/dv`, `e = tau_1(1) + v tau_0(H)` (S01 C3, MEC s.6). Source: UR s.6.4; TH s.1.3.

**D21 Vanishing covector (hatted).** `f(x) := Psi0^(-1) Yhat^mod(x)`, `Yhat^mod := Lambda' Y^C(v(u), s^4 x) = u tau^(-1/2) y(tau x)`,
`y_b(t) := Y^C_b(1,t)`; `phi_(i,m) := (-1)^m [x^(-m)] f_i(x)` = `(u tau^(-1/2)/5) tau^(-m) zeta^(i(m-2))/(4^q q! (m-q)!)`, `q = floor(m/5)`.
Source: GH s.0, s.1-s.2.
