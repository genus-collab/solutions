# DIGEST_A, part 4: s.2 Statements (continued: remarks A18-A27) and s.3 Proofs

Abbreviations: `DIGEST_A.md` header. Statements A1-A17: `DIGEST_A_3.md`.

## s.2.5 Remarks for the introduction, the natural-form section and the evidence appendix

**A18. `rmk:int-hodge`** (rmk; data to the appendix). Linear Hodge frames fail. At genus 2 every coloured stable graph with an edge
is singular at `q = 5^-5`, in every frame tried (F0 true class, F1 linear Hodge with common weight, F2 exponent +1, F3 weight 5, F4
colour-dependent `Lambda^vee(lambda_i)`, F5 point frame), with the maximal pole `s^-10 = t_c^-2`, also for the three all-genus-0
graphs whose bound `Sum_v (2g_v - 2)_+` is 0 (FQ2 s.3). Reason (exact at genus 2): every decoration has valuation exactly `-10(g-1)`;
the pole comes from the row growth `psi_n ~ L^(4n)`, not from `Delta`; a frame change by a q-independent scalar shifts a row coefficient
only at relative order `s^(4k)`, so the coefficients of `s^-10 .. s^-7` are frame independent and equal the point-frame values, non-zero
for all 115 coloured graphs (FQ2 s.3). Further negative results, all genus 2: no constant covariant scalar frame (FQ2 s.6), no
constant diagonal frame (s.7), no constant symplectic colour-mixing frame (s.8: every u-series coefficient of every colour-summed graph
is independent of a constant `M`) regularises the edge graphs; the tautological q-dependent factorisation `M_0` kills the leading
poles but leaves simple poles in the one-edge graphs (FQ2 s.9-s.10; NFS s.4); the covariant-scalar O(u) correction
`(b_1, b_3) = (113/75, -21/125)` regularises genus 2 but has no common solution with genus 3 (FQ2 s.10-s.11). Status: VERIFIED.

**A19. `rmk:nat-notforced`** (rmk). Through genus 3 the graph sums alone do not single out the natural form: besides it there are the
kinematic plane Pi, the invisible orders 4-6 and 17 isolated rays (NF s.2, quoting OU3 s.6d). The natural form is selected by the
reduction principle (uniqueness of `tau`, regular gauge), not by the graph sum (NF s.2; NFS s.0 (iii), s.4: "'Natural' must not be
read as 'forced by the graph sum'"). Status: VERIFIED (genus 3 classification).

**A20. `rmk:nat-analogy`** (rmk, introduction). The factorisation has the shape of the Hodge one: R = (gauge regular at the singular
point) x (R-matrix of the universal local model in its natural variable). Resolved conifold: model `Int exp(-(x - w log x)/z) dx/x`
(Gamma function, Stirling series, natural variable `z/w`, constant in `q`); formal quintic: model `Int exp(-(x^4/4 + v/x)/z) x^m dx`,
natural variable `z` over the canonical coordinate measured from the conifold point, which degenerates there and so depends on `u`
(NF s.2; UR s.7 table). The counterpart of the Hodge factor is `M_0(tau zhat)` together with the regular gauge, not `M_0` alone
(NF s.2; NFS s.4). Status: analogy (context).

**A21. `rmk:nat-G1`** (rmk, optional). `G_1 = zhat X_1`, `X_1 = -(8/225 - w_1/5)[1 (x) H^2 + H^2 (x) 1] - (1/75)[H^3 (x) H^4 + H^4 (x) H^3]`
(Fourier normalisation of the genus-3 code; the second coefficient frame independent), no `H (x) H`; `H (x) H` enters from `u^2` on,
only in the rank-one block `a zhat (H - zhat) (x) (H + zhat)` with `a = -127/15750 - (2/5) beta` (NF s.1 item 4; NFS items 11-13; the
letter `alpha` of the sources for `a` clashes with the frame parameter, s.4 C21). Status: `G_1` PROVED by hand from the G-equation
modulo the hatted R-equation, VERIFIED in 7 frames (NFS item 12); (R1)-(R3) PROVED-mod(Theorem U) (NFS item 11).

**A22. `rmk:nat-ED`** (rmk, optional). The Euler-dilaton element has zero first-order response on one-edge graphs with a vertex of
genus `>= 2`, but is not an exact symmetry of the model theory (genus-1 anomaly `-1/8`) (NFS item 15; NF s.1 item 6, s.3). Status:
VERIFIED genus 2, 3; "general statement conditional on the primary sector at h >= 3" (NFS item 15). Regularity is not a property of
the support rule alone: single K-regular dressings are visible at `u^2` and cancel in the Euler-dilaton combination of `G_2`
(NF s.1 item 6, VERIFIED).

**A23. `rmk:set-genus1`** (rmk, optional). Genus one from the model alone: `dF_1 = -(1/12) du/u` (NFS item 16, VERIFIED);
`[Q^0] d_T F_1 = -3/8` (FQ2 s.2). Sign conventions `<H>_1 = -1/12`, `<e>_1 = +1/8` versus `+1/12`, `-1/8`: S01 C3 (s.4 C10).

**A24. `rmk:set-formal`** (rmk, introduction: scope). The formal quintic is a different theory from the quintic (constant maps at
`g = 2`: formal `F_2(0) = N_{2,0} = -19/9600`, true `-5/144`) (FQ2 s.2). Nothing is claimed about the true quintic or the HKQ
constant (FQ2 s.4; NF s.3 "the non-formal quintic ... not touched").

**A25. `rmk:nat-vertexterm`** (rmk, section on the vertex graph). The vertex term is not the bare model potential: the leafless
vertex is `c_g (u tau^(-3/2) / xhat)^(2-2g)`, and the dilaton leaves of `R'` convert the model coordinate into `t_c` (NF s.2). The
single-leaf bookkeeping `-(g-1) c_g/15` at order `u`, "observed, not derived" in NF s.2 and NFS s.3 item 7, is derived (S01 E7: PROVED
given ML(g), theory Cor. L) and frame independent (S01 C6: `2 gamma_t + 2 w_1 - 4/3 = 1/15`). The vertex law reads `c_g (T/lambda)^(2-2g)`
in the system's own coordinate `T = Int b_0`; `t_c` enters only as `Pi_0` of the leaf point (Lemma T) (S01 s.5). Proof: other digests.

**A26. `rmk:evi-Mindep`** (rmk, evidence appendix). The total is M-independent; `[u^0]` of the vertex graph equals `[u^0]` of the
point-frame total for every `M` with the frozen symbol at `u^0` (since `E_out = O(u)`); so the value `c_g`, and the gap check of the
total, are not evidence for the factorisation, and "the vertex graph carries the polar part" and "all edge graphs are regular" are one
check, not two (NF s.6, NF s.5 genus-4 review). Only the regularity of the edge graphs tests the factorisation (NF s.6).

**A27. `rmk:nat-univ`** (rmk, plan section 9). Universality: canonical gauges of perturbed cyclic systems (not from any
Picard-Fuchs equation) make every edge graph regular, with the vertex graph `c_g (T/lambda)^(2-2g)`; annihilating non-gauge dressings
are regular; violating annihilation or flat regularity breaks regularity (MEC s.2; S01 E12, s.5; VERIFIED genus 2-4). Lemma C: a
change of `B = prod b_a` is a change of coordinate `u = phi(w)`, so B-universality is coordinate covariance (MS2 s.2, PROVED).
B-perturbed systems: VERIFIED genus 2-4 (MS2 s.3 (c)).

## s.3 Proofs

Marks: [W] written out in sentences in the source named; [S] only sketched or stated with a one-line reason; [V] only VERIFIED;
[M] missing in the sources read; [O] owned by another digest (not read here beyond the entry notes).

**A2** (genus-0 data). 1. The I-function from the definition of D1: stated (FQ2 s.1); "follows from the definition in arXiv:1702.06096
eq. (5)" (stage4/evid/referee/REFEREE.md l.52, grep only) [S]. 2. I-function = J-function after the mirror map: a genus-0 mirror theorem,
not in (Q) [M; s.7 item 2]. 3. The fifth-order operator and PF4 from the coefficient recursion of `I`: stated (FQ2 s.1, doc:fq_ring)
[S; elementary, to be written]. 4. `Delta_i`, `kappa`: FQ2 s.1 [V], `q^12` for `kappa_3, kappa_4` (FQ2 s.5); the cyclic form PF4 =
`(1/varpi_0) D (1/Dt) D (u varpi_0^2 (Dt)^2) D (1/Dt) D (1/varpi_0)` is (CYC) (MS5 s.3.4; stage5/critic/CYC_PROOF.md) [O].

**A3** (Frobenius basis). 1. Indicial polynomial `n (n-1)^2 (n-2)` from `p_4 = u`, `p_3(0) = -2`, ..: doc:fq_conifold "checked" [S].
2. Analytic solutions by recursion [S]. 3. `gamma = 7/10` from the `u^-1` obstruction of `g log u + h`: FQ2 s.1 [V]; derived by hand
in the (CYC) proof (MS5 s.3.4) [O]. 4. Isotropy `B(f, g) = 0`: MS5 s.3.4 (`log u -> log u + 1`; the concomitant is a constant,
non-degenerate, alternating form; Pfaffian 2 derived by hand) [O].

**A4** (WKB polynomiality). 1. The hatted equation (HPF) for `Phi_i(zhat, u)`: doc:ref_rows (explicit form) [S]. 2. Its `zhat^(n+1)`
coefficient is `zeta^-i (1-u)(4n - 5 theta_u) ph_n + (terms in ph_0..ph_(n-1)) = 0`, so `ph_n` is a polynomial, unique unless `5 | n`,
where the kernel `u^(4n/5)` of `4n - 5 theta_u` gives one free coefficient and solvability needs a vanishing residue (doc:ref_rows)
[S]. 3. The bound `deg ph_n <= floor(4n/5)` "with trivial resonances": PROVED twice, elementary, cf. arXiv:0710.0889 Thm 4(ii)
(MS4 s.3; the files are Stage-4 worker notes, not read) [O]. 4. `ph_n(0) = F_(0,n)` in every normalisation (MS4 s.3, PROVED three
times) [O].

**A5** ((QRR)). Quoted (MS5 s.2). Sign check, digest arithmetic: `Sum_{a in w_0} a^-1 = Sum_{j=1..4} (1 - zeta^j)^-1 - 1/5 = 2 - 1/5 =
9/5`, so `h_1 = -(B_2/2)(9/5) = -3/20`, against `-11/60` for the weight `+5`; consistent with MS5 s.2 and with the identity row of FQ2 s.9.
Vertex-level derivation: stage5/lit (LIT s.1 per MS5 s.4.3) [O].

**A6** (graph sum = potential). Quoted with A5; calibration FQ2 s.2 [V]. Nothing to prove in the note beyond stating the quote
(s.7 item 3).

**A7** (valuation count). 1. The s-powers of the hatted factors (edge `s^(4(a+a'+1))`, leaf `s^(4d-3)`, vertex
`s^(4 D_v - (2g_v-2+n_v+m_v))`): OU2 s.1, from `Hhat^a = s^(1-c_a) H^a`, `z = s^4 zhat` and the CohFT rescaling
`Omega_{g,n} -> s^-(2g-2+n) Omega_{g,n}` ("referee's derivation") [S: the writers must derive each power from D12; the leaf power uses
`1hat = s 1`]. 2. Dimension count `Sum_v (3g_v - 3 + n_v + m_v) = 3g - 3 - E + m = Sum_e (a_h + a_h') + Sum_l d_l + Sum_v D_v` [W, OU2 s.1].
3. Arithmetic `4 Sum D_v - (2g-2+m) + 4(3g-3-E+m-Sum D_v) + 4E + 4m - 3m = 10(g-1)` [W, OU2 s.1; re-checked here]. 4. With (H) the
rows and `xhat` are u-series, so `G_Gamma` is a u-series [S, OU2 s.1; the code checks it at run time]. 5. Unhatted variant: the
valuations of FQ2 s.3 [W at genus 2; all genera given A4, FQ2 s.3 last sentences].

**A8** (criterion). 1. `T_Gamma = u^(2-2g) G_Gamma` (A7). 2. Reversion of `t_c = u(1 + c_1 u + ..)`, formulas OU2 s.1, OU3 s.1 [W].
3. Triangularity of `dT_k/da_j`: NF s.6 [W, genus 5; the same argument for every g].

**A9** (outer valuation). 1. `Psi` orthogonal drops out of the edges; outer edges `(1 - G(-psi)^T G(-psi'))/(psi + psi')`, leaves
`z (1 - G^-1 1)` (UR s.6.1; OU3 s.3) [S]. 2. `G = 1 + O(u)` makes both `O(u)` (UR s.6.1; NF s.6: `[u^0] E_R = Eh_0`, `tau(0) = 1`) [W].
3. With A7 a decoration is `u^(2-2g+Sum n_e + Sum n_l)` times a series; hence `O(u^(E+2-2g))` and regularity for `E >= 2g-2` [W, UR s.6.1].
4. Only `G_1..G_(2g-3)` enter; with UR's Theorem D (degree bound, the bound only, NFS s.4) the zhat-orders are at most
`floor((10g-13)/4)` [W, UR s.6.1].

**A10** (factorisation). 1. Theorem U applied to the quintic's cyclic system (needs (CYC)) and the model's at `v(u)`: UR s.3 [O].
2. From the gauge to `R~ = G M_0(tau zhat)`, `G = Psi(u)^-1 Ghat Psi(0)`: UR s.4, Corollary F, "PROVED given (H)" [O]. 3. `tau` forced
by equality of the canonical one-forms: UR s.4, NF s.1 item 4 [O]. 4. Normalisation `P_11(0, z) = 1` = (N): final from `D_0 = 1`
(MS5) [O]. 5. Degree bound: UR s.5 (grading, residue relations, degrees), rule W equivalence NFS item 10 and s.5 (bookkeeping, unrefereed)
[O; S for the equivalence].

**A11** (symbol = model). UR s.2; LG s.2-s.4 (LG s.4: "Why it is true: the scaling limit of the Picard-Fuchs operator (PROVED,
unrefereed)", header only). Inputs: (a) = A4; (b) "the conifold-frame map respects D on the row subring, a standard stationary-phase
lemma" (NFS item 4) [M: no citation in the sources read; s.7 item 11].

**A12** (uniqueness of `tau`). (a) [M]: "last step to be written out" (NFS s.3 item 5); (b) [V] (NF s.1 item 2, s.6). Keep as remark.

**A13** (Theorem U). UR s.3 [O]; audit NFS item 2; range check "independent construction regular through `z^80`" (MS5 s.2 table) [V].

**A14** (main theorem). The chain as the entry notes give it; every step [O] (proofs in the other digests).
1. Formal stationary-phase lemmas: changes of variables, Fubini with the factor `beta = (det H_vv(w)/det H_vv(w_0))^(-1/2)`, parameters,
   the exact delta lemma, Stirling as the stationary phase of Euler's integral with `x^m` (L3), homogeneity Lemma H (MS5 s.3.2;
   stage5/lemmas/FORMAL_SP_LEMMAS s.1-s.9, chart/CHART_PROOF_2 s.2; formal proof of (L3) in FORMAL_SP_LEMMAS s.9.3).
2. Theorem Gamma by Theorem LMI, two-variable chart (stage5/chart/referee_adversary/REFEREE.md s.3); the 7-variable chart
   (CHART_PROOF_3 s.3) as a remark (MS5 s.5 item 2).
3. Flat to hatted: `c_n(L) = L^(4n) ph^SP_n(L^-5)`, so `Phi^SP(w; 1) = h(w)` (MS5 s.2 item 3).
4. Theorem P(b): `Phi^SP = Phi^can` (stage4/poly s.6 with the poly referee's s.3; re-derived stage5/pb/PB_AUDIT s.2, new steps
   s.2.6-s.2.9); its hypotheses include (CYC), (S2), the relabelling, L3(iv), Lemma T(b), Lemma H (MS5 s.3.6, O10).
5. The glue `D_0(z) = Phi^can(z, 1)/h(z)` (pb s.3; uses UR s.4, formal WKB uniqueness and (QRR)); hence `D_0 = c(z; 1)/h(z) = 1` and
   `D_i(z) = D_0(zeta^-i z) = 1`.
6. (F_g) for all `g`; then ST2 s.2 steps 1-8: vanishing solutions (Lemmas V, Q, 1); (G); Prop. S; transport (Lemmas E, T, F, Cor. E');
   quintic Theorem R; Theorem S15; ML_K(h) for all h (M1, LD, D(ii), D(iii), M6-M9: Theorem V, Gap Theorem); conclusion.
Where (QRR) enters: only in step 5, in the form `A^QRR_n(q = 0) = h_n` (MS5 s.2); (Q) enters in steps 6 (arXiv:1303.1043 in transport
and S15; arXiv:1509.06954, arXiv:2212.00320, arXiv:math-ph/0702045, optionally arXiv:2502.18115 in ML_K).

**A15** (gap of the total). 1. A14 for every outer graph [O]. 2. The total of the outer sum equals the genus-g potential in the
conifold frame: A6 (quoted), the frame map is a ring map (D11), and the total does not depend on the factorisation (NF s.6, OU3 s.1
composition law) [S]. 3. Sum: polar part `c_g t_c^(2-2g)` [M as a written corollary: `\draftnote`].

**A16-A17, A18-A27.** Remarks; no proof in the note beyond the elementary identities (A16: `c_g` values, digest arithmetic for
`g = 2..6`: `B_4 = -1/30`, `B_6 = 1/42`, `B_8 = -1/30`, `B_10 = 5/66`, `B_12 = -691/2730` give `1/1200, 1/25200, 1/180000, 1/660000,
691/1023750000`). Their numbers go to the evidence appendix (s.5).
