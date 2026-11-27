# DIGEST_A (key read_A): background and the statement

Design-phase digest for the note in `allgenus/note/`. Subject: the formal quintic, its R-matrix and Givental-Teleman graph sum,
the conifold frames, the hatted theory, the natural form `M = M_0(tau(u) zhat)`, the exact statement of Conjecture N(g), the
evidence, the literature named in these sources. No LaTeX. Author-facing only.
Pieces: `DIGEST_A.md` (s.0, s.1 D1-D10), `DIGEST_A_2.md` (s.1 D11-D20), `DIGEST_A_3.md` (s.2 A1-A17),
`DIGEST_A_4.md` (s.2 A18-A27, s.3), `DIGEST_A_5.md` (s.4, s.5 V1-V5), `DIGEST_A_6.md` (s.5 V6-V13, s.6), `DIGEST_A_7.md` (s.7).

**Abbreviations** (paths relative to HVQ = `formalquintic/hodge_vertex_quintic/`).
NF = `NATURAL_FORM.md`; FQ2 = `FORMAL_QUINTIC_genus2.md`; NFS = `natural_form/SYNTHESIS.md`; MEC = `allgenus/MECHANISM.md`;
S01 = `allgenus/SYNTHESIS_STAGE01.md`; MS2, MS4, MS5 = `allgenus/MECHANISM_STAGE2.md`, `_STAGE4.md`, `_STAGE5.md`;
ST2 = `allgenus/stage2/STAGE2_SYNTHESIS.md`; ST2b = `allgenus/stage2/STAGE2_SYNTHESIS_2.md`; OU2 = `OU_EXTENSION_genus2.md`;
OU3 = `OU_EXTENSION_genus3.md`; UR = `natural_form/theory/UNIFORM_REDUCTION.md`; LG = `natural_form/lg/LG_MODEL.md`;
doc:X = the docstring of module X (module docstring; "doc:X.f" = docstring of function or class f).
**Status labels** are those of the sources: VERIFIED (exact finite computation; range stated), ARGUED (NF only: argument written and
audited by a second agent), PROVED, PROVED-mod(X) (NFS), QUOTED, SKETCH, CONJECTURAL, REFUTED. "Final status" means the status after
MS5; MS5 s.2 supersedes MS4 s.2, which updates MS2 s.1, which supersedes MEC s.1, s.3, s.4 (MS2 header, MS4 header, MS5 header).

## s.0 Sources read, and what was not read

**Read in full.**
* NF s.0-s.6, s.6a, s.4 (file table).
* FQ2 s.1-s.11.
* NFS s.0-s.6.
* MEC s.0-s.6; S01 s.0-s.7.
* The entry notes MS5 (s.1-s.5), MS4 (s.1-s.5), MS2 (s.1-s.6), read for the final form of the statement and its inputs.
* Docstrings only (module and function docstrings, no code): `fq_rows.py`, `fq_conifold.py`, `fq_hodge.py`, `fq_graphs.py`,
  `fq_ring.py`, `natural_form/referee/ref_rows.py`, `allgenus/model/ms_rows.py` (all seven exist).

**Read in part** (for conventions and section numbers only).
* OU2 s.1 (the hatted theory and the valuation count) in full; OU3 s.1-s.2 (hatted theory at genus 3, vertex CohFT, the index bug);
  OU3 lines in s.3 and s.6d found by searching for `Psi`.
* UR: section headers; s.6.1 and the start of s.6.2; the lines of s.0, s.4 (Corollary F), s.6.5 and s.7 found by searching for `Psi`.
* ST2 s.2 (the assembled argument, steps 1-8, the LD box, the assembled theorem); ST2b s.6.2-s.6.5 (hypotheses of the quoted
  theorems, the write-up plan, errata).
* Searched only: `allgenus/stage5/lit/LIT*.md`, `allgenus/stage4/evid/referee/REFEREE.md`, `allgenus/stage4/gamma/GAMMA.md`,
  `allgenus/stage5/chart/CHART_PROOF_3.md` for the words "mirror" and "I-function" (one-line hits quoted in s.6, s.7).
* The root bibliography `refGW/refs.bib` (read-only), scanned for the BibTeX keys of the papers in s.6.
* The template note `resolvedconifold/hodge_vertex_conifold/note/`: OUTLINE.md s.0, the header of NOTATION.md and of LABELS.md (for the
  label style only).

**Not read.** LG, UR (beyond the above), `natural_form/referee/REFEREE_NAT.md`, `CONIFOLD_SYMBOL_genus2.md`, the scan notes
(SCALAR, DIAG, SYMPLECTIC), OU2 s.2-s.8, OU3 s.3-s.9 (except the lines named), the genus-4 and genus-5 verification notes
(`natural_form/genus4_referee/`, `natural_form/genus5_cluster/*`; known here only through NF s.5-s.6), every worker and referee note of
`allgenus/` (gauge, univ, model, tr, xyswap, theory, stage2/4/5 workers and referees), the stage syntheses beyond the parts named, all
code beyond docstrings, all logs, and every cited paper. Everything about Theorem U, (G), transport, Theorem R, the Model Lemma, S15,
(CYC), Theorem P(b), Theorem Gamma and the glue is recorded here only at the level of the entry notes, as dependencies of the
statement; its proofs belong to the other digests.

## s.1 Definitions (part 1: D1-D10)

Format: **D-number. Proposed name.** Definition in final form. *Notation.* *Source.*

**D1. The formal quintic.** The torus `T = (C^*)^5` acts on `P^4` with weights `lambda_0, .., lambda_4` on `H^0(P^4, O(1))`. The
formal quintic invariants are `N_{g,d} = Int_{[Mbar_g(P^4, d)]^vir} e(R pi_* f^* O(5))`, with the equivariant Euler class in the
numerator, defined after (virtual) localisation, and then specialised at `lambda_i = zeta^i lambda`, `zeta^5 = 1` primitive;
`lambda = 1` by degree-0 homogeneity. No `C^2` factor (unrefined). At the fixed point `p_i`: `H|_{p_i} = lambda_i`,
`O(5)|_{p_i} = 5 lambda_i`, tangent weights `lambda_i - lambda_j` (`j != i`).
The genus-g potential `F_g` is the generating series of the `N_{g,d}` in the mirror variable (`F_2(0) = N_{2,0}` in FQ2 s.2, so `d = 0`
is included). *Notation:* `N_{g,d}`, `F_g`, `zeta`, `p_i`, `lambda_i`. *Source:* FQ2 s.1 ("Formal quintic of Lho-Pandharipande
(arXiv:1702.06096 s.0.2)"); MS4 s.2 cites arXiv:1702.06096 eq. (5) for the definition; doc:fq_hodge (`wt`, `W0`). The exact form of the
generating series (variable `Q = e^T`, constant maps) is not written in the sources read: s.7 item 1.

**D2. Genus-0 data** ("all verified exactly", FQ2 s.1).
* I-function at `p_i`: `I_i(q, z) = z q^{lambda_i/z} Sum_{d >= 0} q^d prod_{k=1}^{5d} (5 lambda_i + k z) / prod_{j=0}^{4} prod_{k=1}^{d}
  (lambda_i - lambda_j + k z)` (FQ2 s.1; same formula at `lambda_i = 1` in CHART_PROOF_3 s.3.7, grep only). At `lambda = 1`,
  `prod_j (lambda_i - lambda_j + kz) = (lambda_i + kz)^5 - 1` (stage4/evid/referee/REFEREE.md l.51-52, grep only).
* `C_0 = Sum_d (5d)!/(d!)^5 q^d`; `I_1 = C_0 log q + Itilde_1`; mirror map `T = I_1^Q / I_0^Q` ("(LP)"); `tau_0 = 0` at the zeta point
  (FQ2 s.1; doc:fq_ring).
* `D = q d/dq`, `L = (1 - 5^5 q)^(-1/5)`, `D L = (L^6 - L)/5`, `q = (1 - L^-5)/3125` (doc:fq_ring); `C_1 = D T`.
* Yukawa coupling `Y_3 = <<H, H, H>> = 5 L^5 / (C_0^2 C_1^3)`; its instanton part reproduces the quintic's genus-0 numbers (FQ2 s.1
  lists 2875, 609250, 317206375, 242467530000, 229305888887625 for `d <= 5`; the wording there, "N_{0,d} = 2875, 609250/8 ...", is
  ambiguous between GW and BPS normalisation: s.7 item 20).
* Fifth-order equation (equivariant, `M = z D`): `[M^5 - 1 - q (3125 M^5 + 9375 z M^4 + 10625 z^2 M^3 + 5625 z^3 M^2 + 1370 z^4 M +
  120 z^5)] I = 0`, with symbol `ell^5 (1 - 5^5 q) = 1`: eigenvalues `ell_i = zeta^i L` of `H *` with respect to `D`, and
  `mu_i = zeta^i L / C_1` with respect to `d/dT` (FQ2 s.1; doc:fq_rows).
* Fourth-order Picard-Fuchs operator PF4 (acting on `C_0`, `I_1`): `(1 - 3125 q) D^4 - 6250 q D^3 - 4375 q D^2 - 1250 q D - 120 q`
  (doc:fq_ring, doc:fq_conifold); formally self-adjoint for `D` (doc:fq_ring.concomitant; MS5 s.3.4).
* `Delta_i = 1/eta(e_i, e_i) = C_0^2 zeta^(3i) / L^2`, which tends to 0 at `q_c = 5^-5` (FQ2 s.1). Product in the flat basis:
  `H * H^a = kappa_a H^(a+1)`, `kappa = (1, L^5/(C_0^2 C_1^3), 1, C_0/C_1, C_0/C_1)` (doc:fq_rows.Rows; FQ2 s.1 in the form
  `S(H^(a+1)) = z d_T S(H^a)/kappa_a`). FQ2 s.5: `kappa_3 = kappa_4 = C_0/C_1` rests on `q^12` evidence there; the all-orders cyclic form
  used by the programme is (CYC) (MS5 s.3.4, PROVED; other digest).
* All five colours degenerate together at `q_c`; the canonical coordinates stay finite and distinct (FQ2 s.1).

**D3. Points and coordinates on the q-line.** `u := 1 - 5^5 q`. Conifold point `q_c = 5^-5` (`u = 0`); large-radius point `q = 0`
(`u = 1`, `L = 1`); orbifold point `q = oo` (MS4 s.3). `s := u^(1/5)`, so `L = 1/s` (doc:fq_conifold: "L -> 1/s with u = s^5");
`D = -(1 - u) d/du = -((1 - s^5)/5) s^-4 d/ds` (doc:fq_conifold; FQ2 s.1). Rescaled Givental variable `zhat := z / s^4 = L^4 z`
(NF header; FQ2 s.9). *Notation clash:* `s` is also used for `dT_P - v` (S01 E5), for `s(z)` (MS4 s.3) and for the Hodge-class
coefficients `s_k` (doc:fq_hodge): s.4 item C21.

**D4. Frobenius data at the zeta point.** Flat basis `1, H, .., H^4`. Pairing `eta(a, b) = Sum_i a_i b_i 5 lambda_i /
prod_{j != i}(lambda_i - lambda_j)`, hence `eta(H^a, H^b) = 5 delta_{a+b = 3 mod 5}` (FQ2 s.1); `eta^{ab} = (1/5) delta_{a+b = 3 mod 5}`,
pairs `(0,3), (1,2), (2,1), (3,0), (4,4)` (doc:fq_rows.Rows). Canonical idempotents `e_i`; normalised canonical basis
`f_i = e~_i = Delta_i^(1/2) e_i`, `eta(f_i, f_j) = delta_ij` (doc:fq_rows.Rows; OU2 s.1). `Delta_i^(1/2) = C_0 zeta^(3i/2) / L` with the
convention `zeta^(1/2) := zeta^3` (doc:fq_rows.Rows), i.e. `Delta_i^(1/2) = zeta^(4i) x`, `x = s f_0(u)` in the frame `(0,0)` (OU2 s.1;
in a general frame `x = s varpi_0^cf`, since OU2 s.1 sets `xhat := x/s = varpi_0^cf(u)`). `W_i^0 = prod_{j != i}(lambda_i - lambda_j) /
(5 lambda_i) = zeta^(3i)` at the zeta point (FQ2 s.1; doc:fq_hodge.W0).

**D5. Canonical coordinates.** `D u_i = zeta^i L` (FQ2 s.1 "Rows"; doc:fq_rows); `u_i^c` = the value at `u = 0`. The natural variable
of NF s.0 is `tau(u) zhat = -(5/4) zeta^i z / (u_i - u_i^c)` (D14). [Digest check, not a source: `D = -(1-u) d/du` and `L = u^(-1/5)`
give `u_i - u_i^c = -zeta^i Int_0^u t^(-1/5)(1 - t)^(-1) dt = -(5/4) zeta^i u^(4/5) Sum_n 4u^n/(5n + 4)`, consistent with D14.]

**D6. The R-matrix and its rows** (doc:fq_rows.Rows; conventions there cite FRAMEB s.0.1-0.2, A2-A3, not read).
* `R_{ai}(z) = e^(-u_i/z) [S(H^a)]_i`, normalised by `R|_{q=0} = diag(h_{wt_i})` in the f-basis, i.e.
  `R_{ai}(q=0) = h_{wt_i}(z) zeta^(ia) (W_i^0)^(-1/2)`.
* `P_a^(i)(z) := Delta_i^(1/2) R_{ai}(z)`; `P_a^(i) = (C_0/L) X_a`, `X_0 = psi_i/C_0`, `X_{a+1} = (1/(C_1 kappa_a)) (z D + zeta^i L) X_a`;
  so `P_0^(i) = psi_i / L` and `P_a^(i)(q=0) = zeta^(ia) h_{wt_i}(z)`.
* Symplecticity: `Sum_{ab} P_a^(i)(z) eta^{ab} P_b^(j)(-z) = delta_ij Delta_i`. Z/5 covariance: `P_a^(i)(z) = zeta^(ia) P_a^(0)(zeta^(-i) z)`
  (FQ2 s.1).
* Normalised rows `Q := P / Delta^(1/2)` (doc:ref_rows: "Q = P/Delta^(1/2) = zeta^i X_a"); OU2 s.1 writes `Q[a][i][n] = [z^n]
  eta(H^a, R(z) e~_i)` with `val_s Q_{a,n} = -4n + c_a - 1`, and FQ2 s.3 gives `val_s P_{a,n}^(i) = -4n + c_a` after the frame map
  (index placement: s.4 item C4).
* WKB form (doc:fq_rows.WKB): `I_i = e^(u_i/z) z psi_i`, `psi_i = Sum_n z^n psi_{i,n}`, `psi_{i,0} = L`; for `n >= 1`,
  `D(psi_n/psi_0) = r_n` is a Laurent polynomial in `L`; the residues of `5 r_n/(L(L^5 - 1))` at `L = 0` and at the fifth roots of unity
  must vanish (no logarithm); `phi_n = psi_n/L` is fixed by its value `phi_n(L = 1) = [z^n] h_{wt_i}(z)` (the q = 0 normalisation).
* The same in `(zhat, u)` (doc:ref_rows): `psi_i = L Phi_i(zhat, u)`, `Phi_i = Sum_n zhat^n ph_n(u)`, `ph_0 = 1`; the equation (HPF)
  makes `ph_n` a polynomial in `u`, unique unless `5 | n`, where the coefficient of `u^(4n/5)` is free and is fixed by
  `ph_n(u = 1) = [z^n] h_{w_i}(z)`; for `5 not| n` the identity `ph_n(1) = h_n` is a check. Rows: `Qh_0 = (zeta^i/C0) Phi_i`,
  `Qh_{a+1} = kh_a A_{e_a} Qh_a`, `e = c - 1 = (-1, -2, 2, 1, 0)`, `kh = (1/C1, C0^2 C1^2, 1/C1, 1/C0)`,
  `A_e = zeta^i + zhat (1 - u) ((4/5) theta_zhat - e/5 - theta_u)`, `C0 = varpi_0` (frame), `C1 = D t_c`.

**D7. QRR weights and the normalising series.** `w_i = wt_i = {lambda_i - lambda_j}_{j != i} U {-5 lambda_i}`; at the zeta point
`w_i = {zeta^i - zeta^j}_{j != i} U {-5 zeta^i}` (doc:fq_hodge.wt: "numerator twist enters with a MINUS sign"; doc:ref_rows; MS5 s.2).
`h_w(z) = exp(-Sum_{k >= 1} B_2k/(2k(2k-1)) z^(2k-1) Sum_{a in w} a^(1-2k))` (doc:fq_hodge.hv, citing FRAMEB 0.1; doc:ref_rows; MS4 s.1),
written `prod_a Theta(-z/a)` in FQ2 s.1 and doc:fq_hodge. `h := h_{w_0}`, `h_n := [z^n] h(z)`; `h_1 = -3/20` (MS5 s.2).
**(QRR)**: `R(q = 0) = diag(h_{w_i})`, used only in the form `A^QRR_n(q = 0) = h_n` (MS5 s.2). Code form: D6 (`phi_n(1) = h_n`,
`ph_n(1) = h_n`, `hcoef(p, zeta, i, nmax) = [z^n] h_{w_i}(z)`, doc:ref_rows).

**D8. Hodge class; localisation vertex.** `Hodge_w := prod_a Lambda^vee(w_a)/w_a^g = exp(Sum_{k odd} s_k ch_k(E))`,
`s_k = -(k-1)! Sum_a w_a^(-k)` (doc:fq_hodge.class_s, citing FRAMEB A3.2; `w = None` gives the point class). Constant-map vertex at
`p_i`: `(W_i^0)^(g-1) Hodge_{wt_i}` (FQ2 s.1). Vertex integrals `Int_{Mbar_{g,n}} Hodge_w prod psi^d` by Mumford's GRR formula in the
Faber-Pandharipande form (doc:fq_hodge.VertexTable, `_hi`); Witten-Kontsevich numbers by string, dilaton, DVV (doc:fq_hodge.wk).
The convention for `Lambda^vee(a)` is not written in the sources read (s.7 item 4).

**D9. Givental-Teleman graph sum (point frame).** Leg-free stable graphs of genus g; a coloured graph has a colour `i` per vertex.
Vertex: the TFT at colour `i` (normalisation `Delta_i`) with psi-classes integrated (Witten-Kontsevich numbers). Edge between colours
`i, j`: `E_ij(z, w) = [delta_ij Delta_i - Sum_{ab} P_a^(i)(-z) eta^{ab} P_b^(j)(-w)] / (z + w)`, coefficients
`E_{mn} = Sum_{k=0}^{n} (-1)^k G_{m+1+k, n-k}`, `G_{mn} = delta_ij Delta_i delta_{m0} delta_{n0} - (-1)^(m+n) Sum_{ab} eta^{ab}
P_{a,m}^(i) P_{b,n}^(j)` (doc:fq_rows.edge_coeffs, citing FRAMEB C Step 1 and "draft (1.1)"). Dilaton (kappa) leaves:
`Dil_i(d) = (-1)^d [z^(d-1)] P_0^(i)`, `d >= 2` (doc:fq_rows.Dil). Decorations: a psi-exponent per half-edge and dilaton-leaf exponents
`d >= 2` per vertex, grouped into Aut-orbits (doc:fq_graphs.Decorated). Hodge frames: rows multiplied by a z-series, e.g.
`R' = R h_{(lam)}^(-1)` with Hodge-class vertices (doc:fq_rows.rows_times; FQ2 s.3, frames F0-F5).
Model form of the same factors (doc:ms_rows, basis `h^a := N_a H^a`, `N = (1, -1, -1, 1, 1)`): TFT `omega_{g,n}(h^b1, .., h^bn) =
5 [Sum b_i = 3 - 3g mod 5]`; legs `R^-1_k h^a = (-1)^k g^a_k h^(a-k)`; leaves `T(e) = -R^-1_{e-1} h^0 = (-1)^e g^0_{e-1} h^(1-e)`
(`e >= 2`, carries `psi^e`); edges `E(k, l) = [psi^k psi'^l] (eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T)/(psi + psi')`. The R-action used
by the proof is the one of arXiv:1303.1043 Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9 (edge formula, composition, translations; ST2 s.2).

**D10. The differential ring N** (doc:fq_ring). Generators `c_j := D^j C_0`, `i_j := D^j Itilde_1` (`j = 0..3`), `Lam := log q`
(`D Lam = 1`; physical elements are Lam-free), `L^(+-1)`; `D^4` reduced by PF4 (`c_4`, `i_4` as in doc:fq_ring); denominators `c_0`
and `W := c_0^2 C_1 = c_0^2 + c_0 i_1 - c_1 i_0`; the Wronskian relation `B(varpi_0, varpi_1) = 0` (Lagrange concomitant of the
self-adjoint PF4) is linear in `i_3` and eliminates it. Rows, edges, leaves and graph terms are elements of N (FQ2 s.1 "Rows";
FQ2 s.4: `F_2` has 4256 terms, denominators `C_0^61 C_1^15`, of the Lho-Pandharipande shape
`C[L^(+-1)][A_2, A_4, A_6, C_0^(+-1), C_1^(-1), K_2]`).
