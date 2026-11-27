# DIGEST_C1: the reduction to the model, first part (piece 1 of 6)

Key `read_C1`. Design-phase digest for the note in `allgenus/note/`; no LaTeX. Scope: the vanishing tower, the slice K, frames;
the transport lemmas (E, T, F, Lemma 0, Cor. E'); edges in `K_0 (x) K_0` and the leaf point in `-z + K`; quintic Theorem R; the
model's singular structure; universality (Cor. U, B-universality via Lemma C) and the failing controls; what TR_MODEL and
XYSWAP_MODEL establish and which part the final chain uses. This digest is author-facing: file names, statuses and adjudications
below must NOT enter the note (binding rule 5); they are here so that the writers can apply the corrections.

Pieces: this file (s.0, s.1: definitions D1-D19); `DIGEST_C1_2.md` (s.2 parts A-C: tower and slice, transport, quintic Theorem R);
`DIGEST_C1_3.md` (s.2 parts MOD, TR, UNI: model singular structure, TR and x-y swap on the model, universality); `DIGEST_C1_4.md` (s.3
proofs, step by step with locations); `DIGEST_C1_5.md` (s.4 conventions CV1-CV18, s.5 evidence EV-*); `DIGEST_C1_6.md` (s.6 literature,
s.7 open issues O-1..O-22).

**What the final chain uses from these sources** (MECH5 s.2 -> S5A S8 -> S2A s.2 steps 4, 5, 7; details in s.2): theory Lemma 0, Lemma E,
Cor. E', Lemma T, Lemma F and Theorem R (a), (b) (RD2 s.3.1-s.3.2; second write-up MS2 s.6-s.7); the coefficient count M5 (S2A);
TRM s.1-s.2 (TR = model CohFT via arXiv:1509.06954 Thm 8.1, `R_TR = M_0`) and TRM s.4.1 (Gap Theorem, route 1); XY s.1 + s.3
(L-regularity via arXiv:2212.00320 Thm 7.1) and XY s.2.3 (Gap Theorem, route 2); the kernel dictionary (S2B K1). Universality (theory
Cor. U, BVRF Lemma C) and the univ/model experiments are not chain steps: s.9 of the note and the evidence appendix.

**Source abbreviations** (paths relative to `HVQ/allgenus/`, HVQ = `formalquintic/hodge_vertex_quintic`):
* RD1, RD2, RD3 = `theory/REDUCTION_DRAFT.md`, `_2.md`, `_3.md`; TRF = `theory/referee/REFEREE.md`.
* MS1, MS2 = `model/MODEL_SINGULAR.md`, `_2.md`; MRF = `model/referee/REFEREE.md`.
* UN1, UN2 = `univ/UNIVERSALITY.md`, `univ/UNIVERSALITY_E3_E5.md`; URF = `univ/referee/REFEREE.md`.
* TRM = `tr/TR_MODEL.md`; TRRF = `tr/referee/REFEREE.md`. XY = `xyswap/XYSWAP_MODEL.md`; XYRF = `xyswap/referee/REFEREE.md`.
* SYN01 = `SYNTHESIS_STAGE01.md` (Stage 0-1 critic: E1-E13, C1-C11, s.7 overclaims); MECH1 = `MECHANISM.md`; MECH2, MECH4, MECH5 =
  `MECHANISM_STAGE2.md`, `_STAGE4.md`, `_STAGE5.md`.
* S2A, S2B = `stage2/STAGE2_SYNTHESIS.md`, `_2.md` (M1-M10, steps 1-8, K1-K11, O1-O27); BVRF = `stage2/bvar/referee/REFEREE.md`;
  S5A = `stage5/STAGE5_SYNTHESIS.md` (chain S1-S9).
* Cited but not read by me: GI = `gauge/GAUGE_IDENTITY.md` (Lemmas V, Q, 1; (G); Prop. S); UR = `natural_form/theory/UNIFORM_REDUCTION.md`
  (Theorem U, (R1)-(R3)); NF = `../NATURAL_FORM.md`; LG = `natural_form/lg/LG_MODEL.md`.

Status labels are those of the sources: PROVED (argument written, or a quoted theorem with arXiv number), VERIFIED (exact computation,
range stated), CONJECTURAL, SKETCH, REFUTED; referee verdicts CONFIRMED / WEAKENED. "Final status" = after Stage 2 (S2A, S2B) and
Stage 5 (MECH5, S5A).

## s.0 Sources read

**Read in full** (all sections): RD1, RD2, RD3, TRF; MS1, MS2, MRF; UN1, UN2, URF; TRM, TRRF; XY, XYRF; SYN01; MECH1, MECH2, MECH4,
MECH5; S2A; S2B.
**Read in part**: BVRF s.0 (corrections 1-3) and s.3 (Lemma C, statement, proof, hatted encoding); S5A s.2 (S2-S9) and s.3 (lines
85-130); a text search of all `stage5/**/*.md` for the items of this digest (Theorem R, Lemmas E, T, F, 0, K_0, ML_K, Cor. U,
Lemma C, universality): Stage 5 does not change any of them; it only cites them as "transport, quintic Theorem R" in S8 (S5A s.2).
For the style: the first 60 lines of `resolvedconifold/hodge_vertex_conifold/note/LABELS.md` (label prefixes by section file).
**Not read**: GI, UR, NF, LG (they supply the inputs (G), Prop. S, Theorem U, the homogeneity of the model and the definition of
Conjecture N; digested by other readers); the Stage-2 worker notes (DICT_PROOF, LOGLEG_PROOF*, FALSIFIER_B, S15_PROOF, BVAR) and
their referees except BVRF s.0, s.3; `gauge/GAUGE_IDENTITY_HAT.md` (cited for "Theorem D", the polynomiality of hatted edges);
all code and logs (only quoted through the notes); the cited papers (their statement numbers are reported as the sources give them,
s.6). The conifold note's OUTLINE, NOTATION, DRAFT_NOTES and main.tex were not read by me.

## s.1 Definitions (final form; each with notation and source)

**D1. State space, pairing, charges.** `H = Q^5` with flat basis `phi_a = H^a`, `a in Z/5` (`phi_0 = 1` the unit, `phi_1 = H` the
residue direction); `eta(phi_a, phi_b) = 5 delta_(a+b = 3 mod 5)`, `eta^-1 = (1/5) Sum_a H^a (x) H^(3-a)`; charges
`c = (c_0, .., c_4) = (0, -1, 3, 2, 1)`. Source: RD1 s.1.1; MS2 s.9 (eta^-1); TRM s.2.2.

**D2. Loop space, symplectic form.** `Hc = H((z^-1))`, coefficients in `Q[[u]]` where needed; `Omega(f, g) = Res_z eta(f(-z), g(z)) dz`;
`Hc_+ = H[[z]]`, `tau_k(phi_a) := z^k phi_a` (an ancestor insertion `psi^k phi_a`); `Hc_- = z^-1 H[[z^-1]]`. The covector `Omega(f, .)`
pairs `tau_m(e)` with `eta([z^(-m-1)] f(-z), e)`. Solutions are written in LOWERED components `Y_a = eta(phi_a, J)`;
`<Y, n> := Sum_a Y_a n^a = eta(J, n)`. Lowered and vector forms are interchangeable (`eta = 5 J`, `J` the involution `a <-> 3 - a`,
`P^vec = J P J`; (Sp) holds in either form). Source: RD1 s.1.1; RD2 s.3.1 (Lemma E); TRF s.2 (Lemma T paragraph).

**D3. The model (vertex theory).** For `v != 0`, `Omega^mod_v` = the Givental-Teleman CohFT of the Landau-Ginzburg model of the
Hurwitz Frobenius manifold `H_(0;1,4)` on its Z/5-symmetric residue line, at the point `v`; equivalently the frozen theory
`M_0.omega^` rescaled by homogeneity (PROVED in LG s.3, UR s.1; not read). Its correlators `<prod_i tau_(k_i)(phi_(a_i))>_h(v)` are
ANCESTOR correlators (integrals over `Mbar_(h,n)` with its own psi classes) and homogeneous: `= C v^e`,
`e = 2 - 2h + Sum_i w(k_i, a_i)`, `w(k, a) = (4k + c_a - 4)/5`, `C = 0` unless `e in Z` (Z/5 selection). **Regular**: `e >= 0`
(equivalently the value lies in `Q[v]`); **singular**: `e < 0` and `C != 0`. The curve is written in three normalisations (s.4 CV8).
Source: RD1 s.1.2; MS1 s.1 (w with the opposite sign, s.4 CV3); S2A s.1 ("lies in Q[v]").

**D4. Ancestor potential, base point, dilaton.** `Fb_h(v; q) = Sum_n <(q + z phi_0)^n>_h / n!` for `q in Hc_+` (dilaton shift; base
point `q = -z phi_0`). Homogeneity: `Fb_h(lambda q) = lambda^(2-2h) Fb_h(q)` (`h != 1`), `Fb_1(lambda q) = Fb_1(q) - (5/24) log lambda`.
Correlator form: `<tau_1(1) X>_h = (2h - 2 + n) <X>_h`. Source: RD1 s.1.2; MS1 s.1.

**D5. Free energies, vertex constants.** `F_h = c_h v^(2-2h)` (`h >= 2`), `c_h = |B_2h| / (2h (2h-2) 5^(h-1))
= (-1)^(h-1) B_2h / (2h (2h-2) 5^(h-1)) = (-5)^(1-h) B_2h / (2h (2h-2))`; `F_1 = -(1/12) log v`; `F_0 = -(5/2) v^2 log v` (only
derivatives of order `>= 3` are compared; normalisation `<HHH>_0 = -5/v`). `Phi_h(x)` denotes the same functions of a variable.
Values `c_2..c_5 = 1/1200, 1/25200, 1/180000, 1/660000`; `c_6 = 691/1023750000` (a prediction in SYN01 C2). Source: RD1 s.1.2, s.2;
MS1 s.0; S2A s.1; SYN01 C2 (four equal forms, checked `g = 2..12`).

**D6. Vanishing solutions.** `Y^C(v; z)`: the pure-`1/z` (residue, vanishing-cycle) solution of the model's system, lowered form;
`Y^C_A(u; z)` the same for a self-dual cyclic one-pole system `A`; `Y^C_Q` for the formal quintic; `J^C := eta^-1 Y^C`. Properties
(PROVED in GI s.2, Lemmas V, Q; quoted by RD1 s.1.3): both lie in `Q[z^-1][[u]]^5`, equal `-e_1` at `u = 0`; `Y^C_b` has only the
powers `z^-n` with `n = c_b + 1 mod 5`; `[z^-1] Y^C_0 = -v`, `[z^-1] Y^C_(A,0) = -T_A`, `T_Q = t_c`. Leading terms
`Y^C = (-v/z, -1, -v^4/(24 z^4), -v^3/(6 z^3), v^2/(2 z^2)) + ..` (TRM s.3). Decay: `[z^-n] Y^C = O(v^(n - floor(n/5)))` (RD2 s.3.1).
Model form: `Y^C_b(z) = -(1/5) eta(S(z) H^2, H^b)`, `S(z) = 1 + S_1/z + ..` the homogeneous calibration (`z d/dv S = (H *) S`,
integration constants fixed by homogeneity); the column `S(z) H^2` is log-free, only `S(z) H` has a logarithm (`S_1 H = -log(v) H^2`).
Source: RD1 s.1.3; MS1 s.1; MRF item 5 ("log-free", not "no resonance").

**D7. The vanishing tower.** For `q in Hc_+`, `i >= 0`: `Pi_i(q) := (1/5) [z^i] eta(S(-z) H^2, q(z))` (= `-[z^i] <Y^C(-z), q(z)>` by D6);
`Pi_0 =: PiC`, the vanishing covector. Shift rule `Pi_i(tau_k(H^b)) = PiC(tau_(k-i)(H^b))` (zero for `k < i`). Meaning: `Pi_i(q)` is
the descendant time `t^H_i` of the ancestor vector `q` (descendant times `= [S^-1(z) q(z)]_+`, `S^-1(z) = S^*(-z)`). Values:
`PiC = 1` on `tau_0(H)`, `-v` on `tau_1(1)`, `-v^2/2` on `tau_2(H^4)`, `-v^3/6` on `tau_3(H^3)`, `v^4/24` on `tau_4(H^2)`,
`v^4/96` on `tau_5(H)`, `-v^5/480` on `tau_6(1)`; `PiC(tau_k(H^b)) = 0` unless `b + k = 1 mod 5`; `PiC(tau_k(H^b)) =
(-1)^(k+1) [z^-k] Y^C_b`; SYN01 C4 records a closed form `1/(4^q q! (m-q)!)` with support rule `b = 1 - m` (indexing in GI Lemma V,
not read). Theory's flag: `V^(k)(q) := (-1)^k [z^(k-1)] <Y^C(-z), q(z)> = Omega(z^-k J^C, q)`, `k >= 1`; by the two definitions
`V^(k) = (-1)^(k+1) Pi_(k-1)` (s.4 CV1). The "D-directions" (D-tower) are `tau_k(H)`, `k >= 1`. Source: MS1 s.1; RD1 s.1.3; SYN01 C4;
S2A M3 (`Pi_0 = PiC := (-1)^(n+1) [z^-n] Y^C_b`).

**D8. The slice K and K_0.** `K := {q in Hc_+ : Pi_i(q) = 0 for all i >= 1}` (all descendant times `t^H_k`, `k >= 1`, vanish);
`K_0 := K cap ker Pi_0 = {n in Hc_+ : <Y^C(-z), n(z)> in z^-1 Q[[z^-1]]}` (theory's `T_v`); `K = K_0 (+) Q tau_0(H)`. Affine slice
`-z + K = {q in -z phi_0 + Hc_+ : Pi_i(q) = 0, i >= 1} = -z phi_0 + (K_0 (+) Q phi_1)` (theory's `Lambda_v`), with
`Pi_0(-z phi_0 + n + s phi_1) = v + s` for `n in K_0`. Coefficients: `Q[v]` for homogeneous elements, `Q[[u]]` after `v = v(u)`.
Notation of record: `K` (the Stage-2 conventions call it `K_vc` to avoid UR's frame `K = (1, uH, H^2, H^3, H^4)`), `K_0`. Source:
MS1 s.1; RD1 s.1.3; S2A preamble; SYN01 C9.

**D9. Frames of K_0.** For `(k, a)`, `a != 1`: `frame(k, a)` := the unique element `z^k phi_a + Sum_(j<k) f_j z^j phi_1` of `K_0`,
`f_j = const . v^((4(k-j) + c_a + 1)/5)` (exponent `>= 1`; non-zero only for integral exponent). At `v = 1`: `frame(1,0) = z + H`
(general `v`: `z phi_0 + v phi_1 = e`), `frame(2,4) = z^2 H^4 + H/2`, `frame(3,3) = z^3 H^3 + H/6`, `frame(4,2) = z^4 H^2 - H/24`;
`frame(k,a) = z^k H^a` for `(0,0), (0,2), (0,3), (0,4), (1,2), (1,3), (1,4), (2,3)`; the first frame that depends on the `z^-5` terms
of `Y^C` is `frame(6,0) = z^6 phi_0 + v z^5 H - (v^5/120) H`. No frame has `a = 1`. Source: RD1 s.1.3; TRF s.0 item 3(ii); S2A M5.

**D10. Model basis of K.** `kappa(X, k) := [S(z) X z^k]_+ = Sum_(j<=k) S_j X z^(k-j)`, `X in {1, H^2, H^3, H^4}` (the ancestor images of
the descendant insertions `tau_k(X)`); coefficients are non-negative powers of `v` (`S_j H^a` carries `v^((c_a - c_(a+j) + 4j)/5)`,
exponent `>= 1` for `j >= 1`, MRF item 8(c)). `K = Q tau_0(H) + span{kappa(X,k)}`, the `kappa` span `K_0`. A basis insertion
`tau_k(H^b)` lies in `K` iff `b = 0, k <= 1`; `b = 1, k = 0`; `b = 2, k <= 4`; `b = 3, k <= 3`; `b = 4, k <= 2`. Source: MS1 s.1.

**D11. Dressings and their hypotheses.** A dressing of the model at `v(u)`: `P(u, z) = Sum_(k>=0) z^k P_k(u)`, `P_k in gl_5(Q[[u]])`
(**flat-regular**), Z/5-**graded** (`P_k = Sum_b p_(k,b) E_(b+k,b)`), `P_0` invertible, acting on lowered components.
**(Sp)** `P(-z)^T eta P(z) = eta`. **(A)** (annihilation of the vanishing period) `Y' := P Y^C(v(u))` has no positive power of `z`.
Put `T_P := -[z^-1] Y'_0` (the transported vanishing coordinate), `d := (P_0)_33` (the factor by which `P_0` rescales the unit in vector
form; the dilaton factor), `R := P_0^-1 P`. Examples: (i) the canonical Theorem-U gauge from the pulled-back model to a self-dual
cyclic one-pole system `A` (`Y' = Y^C_A`, `T_P = T_A`); (ii) `D P_A` with `D` flat-regular, graded, symplectic and
`(D Y^C_A)_(>0) = 0` (non-gauge dressings). `v(u)` is the model coordinate pulled back (D15); `v(0) = 0`, `v'(0) != 0`. Source: RD1 s.1.4.

**D12. Edge, leaf point, graph sum.** For `R(z) = 1 + O(z)` symplectic: `E_R(psi, psi') := (eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T)/(psi + psi')`
(a power series by (Sp)); `E_P := E_R` for `P = P_0 R`. Leaf point `q_T := P^-1(-z phi_0)` (vector form; contains the dilaton part
of the translation `z(1 - R^-1(z) 1)`), `qt := d q_T`. `Omega^P := P.Omega^mod_(v(u))`: the R-matrix action of arXiv:1303.1043 s.2
(edges as above, legs `R^-1(psi)`, translation `z(1 - R^-1(z) 1)`) of `R`, preceded by the constant frame change `P_0`. For a
leg-free stable graph `Gamma` of genus `g`, `Cont_Gamma` contracts the bivectors `E_P` along the edges with the vertex tensors
`d^(n_v) Fb_(h_v)(q_T)` (the `n_v`-th derivative at `q_T`); `F^P_g = Sum_Gamma |Aut Gamma|^-1 Cont_Gamma`. **Vertex graph** = one vertex,
no edge; **edge graph** = a stable graph with at least one edge; **regular** (flat frame) = `Cont_Gamma in Q[[u]]`. Hatted form (not in
this digest): the outer sum `R'.(M.omega^)`, `M = M_0(tau(u) zhat)`; a genus-g term `u^(2-2g) G(u)` is regular iff
`[u^0..u^(2g-3)] G = 0` (UN1 header); flat = hatted graph by graph is Theorem S15 (S2A step 6). Source: RD1 s.1.5; RD2 s.3.1-s.3.2.

**D13. Canonical and stored gauge.** `P_can`: the Theorem-U gauge (pulled-back model to the quintic) normalised by `P_11(0, z) = 1`;
by (G) equivalently `(P Y^C)_+ = 0`. Theorem U leaves `P -> P c(z^5)` free. Stored: `P_stored = P_can (1 + m5 z^5 + ..)`,
`m5 = -43397/807975000` (`P^stored_11(0, z) = 1 - (43397/807975000) z^5 - ..`; `z^10` coefficient `-765368951400317/1305647201250000000`).
Source: MECH1 preamble, s.6; UN1 s.1, s.2.1; URF s.2.1; SYN01 C11.

**D14. The Model Lemma ML_K (final form).** ML_K(h), `h >= 0`: (i) for `n >= 1`, `2h - 2 + n > 0`, `X_1..X_n in K` homogeneous with
coefficients in `Q[v]`: `<X_1..X_n>_h(v) - F_h^(n)(v) prod_i Pi_0(X_i) in Q[v]`; (ii) for `h >= 2`: `F_h = c_h v^(2-2h)` (the gap).
Equivalent forms: (ML1) `<H^n>_h = F_h^(n)(v)`, (ML2) a correlator of K-generators (`tau_0(H)`, `kappa(X,k)`) with at least one
`kappa` has `e >= 0` or vanishes (MS1 s.5); theory's ML_Lambda(h): for homogeneous `e_1..e_n in K_0`, `m >= 0` (`n + m >= 3` if
`h = 0`, `>= 1` if `h = 1`), `<e_1..e_n tau_0(H)^m>_h - delta_(n,0) Phi_h^(m)(v)` is a regular monomial, zero if its exponent is
negative; function form `Fb_h | (-z + K) = Phi_h o Pi_0 + R_h`, `R_h` regular (RD1 s.2). Source: S2A s.1 (final); MS1 s.5; RD1 s.2.
(ML_K is proved in the note's model sections, from S2A s.1 M1-M10; other readers.)

**D15. Self-dual cyclic one-pole systems.** `th Y_a = b_a Y_(a+1)`, `b = (f, B/(f^2 g^2), f, g, g)`, `B = prod b_a`, `ord_0 B = -1`
(`th` as in UR s.1, not read); quintic `B = -1/(u (1-u)^5)`; univ family `f = t_c' F(u)`, `g = -(varpi_0/(1-u)) Gm(u)`. Model
pull-back coordinate `v_X`: `v_X'^5 = -B v_X`, `v_X(0) = 0`, `a := v_X'(0)`, `a^4 = -res_0 B` (quintic `v_q'(0) = 1`). Vanishing
coordinate `T = Int_0^u f`. Residues `r`, `rt` of the two poles; `lambda` by `r = -1/lambda^2`, i.e. `lambda = f(0) |g(0)|` (SYN01 C10).
Source: UN1 s.0-s.1; BVRF s.3; RD2 s.3.3; SYN01 C10.

**D16. Legs, leg class, pairing (model as TR).** A leg is a function `Phi` on the curve; `<omega_(g,n), Phi_1..Phi_n> := Sum
Res_(a_1) .. Res_(a_n) prod_j Phi_j(z_j) omega_(g,n)` (residues at the branch points). Functions of `x` pair to zero (linear loop
equation). Leg class `L = Q[params][z] + Q[params][x]` (XY), on the line `Q[v][p] + Q[v][x]` (S2A M4, M7). Source: XY conventions,
s.0 item 4, s.3; S2A s.1.

**D17. TR data of the model curve.** `x = kappa W`, `W = p^4/4 + u/p` (here `u` is the RESIDUE, not the conifold coordinate; S2A
writes `W_v = p^4/4 + v/p`), `y = p`, `B = dp1 dp2 / (p1 - p2)^2`; final normalisation `kappa = 1`. Branch points `c_j = zeta^-j`
(`u = 1`), `x(c_j) = (5/4) kappa zeta^j`; `Lambda = -5 kappa`; xi-functions `xif_0..xif_4 = (-p^4, 1, p, -p^2, -p^3)/(p^5 - 1)`
(the dual vectors `H^(3-a)/5`), `xif_(a,d) = (d/dx)^d xif_a`; general `v`: `xif_a = eps_a v^(alpha_a) p^(m_a)/(p^5 - v)` (S2A step 7.2).
`F_g^TR := (1/(2-2g)) Sum Res Phi omega_(g,1)`, `dPhi = y dx`. Source: TRM s.2.1-s.2.2, s.4.1; S2A step 7.2.

**D18. Euler-dilaton vector** `e := tau_1(1) + v tau_0(H) = frame(1, 0) in K_0`. Source: RD1 s.1.3; SYN01 C3.

**D19. Isotropic c = 1 sector** `Sigma_C := Q((z^-1)) J^C`; the `Pi_i` are the coordinates dual to its negative half. Source: RD1 s.2.2.
