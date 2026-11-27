# DIGEST_A, part 3: s.2 Statements

Abbreviations: `DIGEST_A.md` header. Proposed label prefixes (by topic): `int` (main statement, introduction), `set` (formal quintic,
genus-0 data, frames, graph sum), `hat` (hatted theory), `nat` (natural form), `mod` (the model), `evi` (evidence appendix). Each entry:
label, type; statement (final form); scope; source; status in the sources; corrections applied; dependencies.

## s.2.1 The formal quintic and its graph sum

**A1. `def:set-fq`** (defn). D1. Source FQ2 s.1. Status: definition, quoted (arXiv:1702.06096 s.0.2 per FQ2 s.1; eq. (5) per MS4 s.2).

**A2. `prop:set-genus0`** (prop). The genus-0 data of D2 at the zeta point: the I-function, PF4 and the fifth-order equation, the
eigenvalues `zeta^i L`, `Delta_i = C_0^2 zeta^(3i)/L^2`, the product `kappa`. Scope: `lambda_i = zeta^i`. Source FQ2 s.1, doc:fq_rows,
doc:fq_ring. Status: "all verified exactly" (FQ2 s.1); `kappa_3 = kappa_4 = C_0/C_1` on `q^12` evidence (FQ2 s.5); the cyclic form of
the quintic system, (CYC), PROVED and refereed twice (MS5 s.3.4). Dependencies: a genus-0 mirror theorem for the formal quintic is
needed to identify the quantum differential equation with PF; no such theorem is in (Q) (s.7 item 2).

**A3. `lem:set-cfbasis`** (lem). PF4 has local exponents `{0, 1, 1, 2}` at `u = 0`; analytic solutions `f_0, f_1, f_2` as in D11; the
logarithmic solution is `g log u + h` with `g = f_1 + (7/10) f_2`; `B(f, g) = 0` for every power-series solution `f`, so every frame
`varpi_0^cf = f_0 + alpha g + beta f_2` is isotropic and `(alpha, beta)` and the scale of `g` are not fixed by the pairing.
Source FQ2 s.1, doc:fq_conifold, MS5 s.3.4. Status: exact computation (FQ2); `gamma = 7/10` and the isotropy PROVED in the (CYC) proof
(MS5 s.3.4: "the automorphism `log u -> log u + 1` of `Q((u))[log u]` makes every power-series solution orthogonal to `g`").
Dependencies: none. Proof owned by the (CYC) digest.

**A4. `lem:set-wkb`** (lem). For each colour `i`, `ph_n(u)` (D6) is a polynomial in `u` of degree at most `floor(4n/5)`; it is
determined by the recursion except, for `5 | n`, its coefficient of `u^(4n/5)`; no logarithm occurs. Equivalently `psi_n/L` is a
Laurent polynomial in `L` of top degree at most `4n`. `ph_n(0)` is the model value in every normalisation. Source doc:ref_rows;
doc:fq_rows.WKB; FQ2 s.1; MS4 s.3. Status: the degree bound "with trivial resonances is PROVED twice (elementary; arXiv:0710.0889
Thm 4(ii)); equality is VERIFIED only" (MS4 s.3); FQ2 s.1 observed "top degree exactly 4n for n <= 7"; "`ph_n(0) = F_(0,n)` in every
normalisation (PROVED three times)" (MS4 s.3). Correction: state the bound, never equality (MS4 s.3).

**A5. `prop:set-qrr`** (quoted input (QRR)). `R(q = 0) = diag(h_{w_i})`, `w_i = {zeta^i - zeta^j}_{j != i} U {-5 zeta^i}`, used only
in the form `A^QRR_n(q = 0) = h_n` (MS5 s.2). In the rows: `ph_n(1) = h_n` (fixes the free coefficient for `5 | n`; a consequence
for `5 not| n`, doc:ref_rows). Status: QUOTED: GJR arXiv:1812.11908v1 eq. (17) (proof of Lemma 5.2) and arXiv:1709.07392v1 eqs.
(19)-(20) at `t = 0`, both justified there by arXiv:math/0008067 Thm 1 with its Remark and arXiv:math/0110142 Thm 1, Cor. 1; a
vertex-level derivation PROVED given its quotes (arXiv:1803.01409v2 eq. (26), Mumford's relation, arXiv:math/0108100v2 Thm 4.1);
uniqueness of R: arXiv:0712.0160v3 Theorem (2) (s.1.9) and Prop. 8.5(ii) (s.8.2), "read only through a summariser" (MS5 s.2).
Twist sign `-5` "settled three ways and pinned by `h_1 = ph_1 = -3/20`" (MS5 s.2). Correction: MS4 s.2 cited arXiv:math/0108100
Thm 9.1, "the untwisted statement"; replaced by GJR (17) (MS5 s.1, s.3.6).

**A6. `thm:set-gt`** (quoted). The genus-g potential of the formal quintic is the Givental-Teleman graph sum of D9 with the actual
(QRR-normalised) R-matrix. Source: implicit in "the actual R-matrix" (MS4 s.2: "It was implicit in Stage 2"); the uniqueness quote of
A5. Status: QUOTED as part of (QRR); calibration VERIFIED against direct localisation (FQ2 s.2). Correction: none in the sources; the
note must state where this enters (s.7 item 3).

## s.2.2 The hatted theory

**A7. `lem:hat-count`** (lem). Let `g >= 2`, `M(zhat)` symplectic. For a decoration of a leg-free genus-g graph with `E` edges, `m`
dilaton leaves, half-edge exponents `a_h`, leaf orders `d_l`, vertex psi-deficits `D_v`, the hatted factors are `s^(4(a+a'+1))` times
the edge factor, `s^(4d-3)` times the leaf factor and `s^(4 D_v - (2g_v - 2 + n_v + m_v))` times the vertex; using
`Sum_v (3g_v - 3 + n_v + m_v) = 3g - 3 - E + m = Sum_e (a_h + a_h') + Sum_l d_l + Sum_v D_v`, the total exponent of `s` is `10(g-1)` for
every decoration. Hence, if the hatted rows and `xhat` are power series in `u` ((H)), `G_Gamma := s^(10(g-1)) T_Gamma` is a power series
in `u` and `T_Gamma = u^(2-2g) G_Gamma`. Scope: every `g`, every symplectic `M` constant in `u` (OU2 s.1); for `M = M_0(tau(u) zhat)` the
same count applies coefficientwise (NF s.5-s.6 use it; digest remark). Unhatted form (FQ2 s.3): after Cf,
`val P_{a,n}^(i) = -4n + c_a`, `val E_{mn} = -4(m+n) - 2`, `val Dil(d) = -4(d-1)`, `val Delta^(g_v-1) = 2(g_v - 1)`, so every decoration
has valuation `-10(g-1)` in `s`, i.e. pole order at most `2g - 2` in `u`. Source OU2 s.1 ("The valuation count"; "re-derived by hand by
the referee (s.7)"), OU3 s.1, FQ2 s.3, ST2 s.2 step 6. Status: argument written; FQ2 s.3 says the all-genus form "needs
`deg_L(psi_n/L) <= 4n` for all n", which is A4 (PROVED). Dependencies: A4, D12, (H) (final: from `D_0 = 1`, D20).

**A8. `lem:hat-regular`** (lem). `T_Gamma` is regular at `u = 0` iff `[u^k] G_Gamma = 0` for `0 <= k <= 2g - 3`. The `t_c`-polar
coefficients are obtained from `a_k = [u^k] G_Gamma` by series reversion; the map is lower triangular with unit diagonal, so
"polar part `= c t_c^(2-2g)`" is `2g - 3` independent affine conditions on the `a_k` given `a_0 = c`. Source OU2 s.1, OU3 s.1, NF s.6
(the Jacobian `dT_k/da_j = [t^k](A^-8 u(t)^j)` at genus 5). Status: elementary (written in NF s.6 for genus 5).

**A9. `lem:nat-val`** (lem). For `R' = Psi(u) G(u, zhat)` with `G = 1 + O(u)`, every outer edge and every outer leaf is `O(u)`. A
genus-g decoration of the outer sum is `u^(2 - 2g + Sum_e n_e + Sum_l n_l)` times a power series (`n_e, n_l >= 1` the u-orders taken from
the factors). Hence an outer graph with `E` edges is `O(u^(E + 2 - 2g))`, every graph with `E >= 2g - 2` is regular, and for the others
only `G_1, .., G_(2g-3)` enter, i.e. by the degree bound `zhat`-orders at most `floor((10g - 13)/4)` of the gauge. Source UR s.6.1;
NF s.6 ("a one-line valuation argument, not evidence"; `[u^0] E_R = Eh_0` and `tau(0) = 1` force `E_out = O(u)`). Status: argument
written (UR s.6.1); "valuation exactly 1" of every entry of `E_out` VERIFIED (1950 entries at genus 5, NF s.6). Dependencies: A7, D15
(`G = 1 + O(u)`; `R' = Psi + O(u)`, OU3 s.3), A10 for the degree bound.

## s.2.3 The natural form and the model

**A10. `thm:nat-fact`** (thm). `R~(zhat; u) = G(u, zhat) M_0(tau(u) zhat)`, `G = 1 + Sum_{n >= 1} u^n G_n(zhat)`, `G_n` a polynomial in
`zhat` of degree at most `floor((5n + 2)/4)`; `M_0` stands on the right. Scope: all orders; every conifold frame. Frame independent:
`M_0`, `tau`, polynomiality, the degree bound, the chain relation; frame dependent: values, sharpness, degrees (NFS item 13). Source
NF s.1 item 1; NFS items 1, 10, 13; UR s.4 Corollary F ("PROVED given (H)"). Status: all orders PROVED-mod((H), Theorem U) (NFS item 1;
NFS s.2 item 2 settles the referee's "CONJECTURAL" in favour of this, Theorem U being the missing existence statement); the degree
bound PROVED-mod((H), canonical normalisation) through the equivalence "rule W = Theorem U + (R1) + (R2) + item 8" (NFS item 10, "index
bookkeeping ... unrefereed", NFS s.5). Final: (H) and the canonical normalisation follow from `D_0 = 1` (MS5 s.1). Corrections: the
degree `floor((5n+2)/4)` is a bound, attained in generic frames (NFS s.2 item 3; NF s.1 item 1: frame `(0,0)` degrees 1, 3, 4, 5, 6, 8,
9, ..); the left placement is false (NF s.1 item 1, VERIFIED). Dependencies: Theorem U (UR s.3), UR s.2 (A11), (H), (N).

**A11. `prop:mod-symbol`** (prop). `M_0` is the R-matrix of the model of D16: `R^_a^(i)(zhat) = N_a zeta^(i(a+1)) J_{c_a}(-4 zeta^-i
zhat)`, with the closed recursion of D16. Source NF s.1 item 3; NFS items 4-5; UR s.2; LG s.2-s.4. Status: VERIFIED (all rows
`zhat <= 18` at two primes, identity row `zhat <= 40`, two stationary-phase engines); all orders PROVED-mod((a) "WKB coefficients are
Laurent polynomials in L", (b) "the conifold-frame map respects D on the row subring, a standard stationary-phase lemma"); independent
of (H) (NFS item 4). Corrections: LG's "Theorem A PROVED" read as PROVED-mod((a), (b)) (NFS s.4); dictionary `k = c_a`, not `x^a dx`
(NFS s.2 item 1). In the final chain the identification enters through ST2 s.2: "(F_g) implies (H) in the window, hence UR s.2's
identification of the conifold symbol with the model's `M_0`". Input (a) is A4.

**A12. `rmk:nat-tau`** (rmk). Uniqueness of the rescaling: (a) at most one `tau` makes `G` polynomial, and it is the canonical
coordinate; (b) in the ansatz `M = M_0(tau zhat)` the regularity of the one-edge graphs fixes `tau` mod `u^(2g-2)`. Source NF s.1
item 2, NF s.6 (controls). Status: (a) ARGUED modulo the hatted R-equation to all orders (NF), "PROVED-mod(hatted R-equation ..., i.e.
Picard-Fuchs + (H))" with "last step to be written out" (NFS item 3, s.3 item 5); (b) VERIFIED genus 2-5. Not used by the main theorem.

**A13. `thm:set-U`** (thm, statement level only). Two cyclic systems `z d/du Y_a = b_a Y_(a+1)` with unit coefficients except for one
simple pole and with equal products `prod b_a` are related by a formal gauge `P(u, z)` whose coefficients are regular at `u = 0`,
unique up to a constant `c(z^5)`; equality of the products is equality of the canonical one-forms and forces `v = u tau^(-5/4)`
(NF s.1 item 4; NFS item 2). Status: PROVED in UR s.3, audited (NFS item 2; MS5 s.2: "(U) ... proved in the project, AI-audited
only"). The exact hypotheses are in UR s.3 (not read); the note proves it (task rule 2).

## s.2.4 The main statement

**A14. `thm:int-main` (Conjecture N(g)).** Let `g >= 2`, fix a conifold frame `(alpha, beta)` and `t_c = g/varpi_0^cf`. Factorise the
hatted R-matrix of the formal quintic as `R = R' M` with `M = M_nat = M_0(tau(u) zhat)` and `R' = Psi(u) G(u, zhat)`. In the outer graph
sum `R'.(M.omega^)` of the genus-g potential: (i) the term of every leg-free stable graph of genus g with at least one edge is regular at
`u = 0` (lies in `Q[[u]]`); (ii) the term of the vertex graph is `c_g t_c^(2-2g) + Q[[u]]`, with
`c_g = (-1)^(g-1) B_2g / (2g(2g-2)) 5^(1-g) = |B_2g| / (2g(2g-2) 5^(g-1))`.
Equivalent forms: `Omega^quintic = P . (Omega^LG at v(u))` with `P` regular (NF s.1 item 7); in the flat frame, every stable graph with
an edge in `P_can . Omega^mod_{v(u)}` lies in `Q[[u]]` and the vertex graph is `c_g t_c^(2-2g) + Q[[u]]` (ST2 s.2 step 5), the hatted and
flat terms being equal graph by graph (Theorem S15, ST2 s.2 step 6).
Scope: per graph (not only the total); `M` exact in `u`; the formal quintic only. Source NF s.1 item 7; ST2 s.2 step 8; MS5 s.2.
Status: NF s.1 item 7: "VERIFIED for g = 2, 3, 4, 5 ...; CONJECTURAL for g >= 6" (superseded); final: holds for every `g >= 2` modulo
exactly (U), (Q), (QRR) (MS5 s.2), (U) proved in the project. Corrections applied in the final form: the s.1.5 identification is not an
input (MS2 s.1: derived as S15 from (F)); (QRR) listed explicitly (MS4 s.2); (F_g) discharged by Theorem Gamma, Theorem P(b) and the glue
(MS5 s.2); (CYC) added to the hypotheses used by P(b) and proved (MS5 s.1, s.3.4, s.3.6); S01 s.7: "Cor. N_9 ... also needs s.1.5
(genus 2 only) and Theorem U" (now covered). Dependencies: (U), (Q), (QRR); the chain MS5 s.2 items 1-6 and ST2 s.2 steps 1-8.
Frame scope: s.7 item 9.

**A15. `cor:int-gap`** (cor, proposed). In every conifold frame the genus-g potential in the conifold frame is
`c_g t_c^(2-2g) + 0 t_c^(3-2g) + .. + 0 t_c^-1 + Q[[u]]`; the constant term and the higher coefficients depend on the frame. Source:
FQ2 s.4 (`g = 2`, four frames), NF s.0, s.5, s.6 (`g = 3, 4, 5`). Status: VERIFIED `g = 2..5`; for all `g` it follows from A14, A6 and
the M-independence of the total (NF s.6), but no source read states it as a corollary: `\draftnote` (s.7 item 10).

**A16. `rmk:int-cg`** (rmk). Values `c_2 = 1/1200`, `c_3 = 1/25200`, `c_4 = 1/180000`, `c_5 = 1/660000`; `c_6 = 691/1023750000`
(prediction). The four forms `(-5)^(1-g) F_g(x = W)`, `|B_2g|/(2g(2g-2)5^(g-1))`, `(-1)^(g-1) B_2g/(..) 5^(1-g)`,
`20^(1-g) (-4)^(3g-3) F^Hurw` agree for `g = 2..12` (S01 C2). `c_g` refers to `g = u + O(u^2)`; the scale of `g` is free
(doc:fq_conifold); `1/1200 = (1/5)(1/240)` and `t_c -> sqrt(5) t_c` would give the usual conifold value; the formal quintic is a
different theory, "so no statement about the HKQ constant is made" (FQ2 s.4).

**A17. `rmk:int-pergraph`** (rmk). The per-graph statement is stronger than any statement about the total: in the stored, non-canonical
normalisation (`P_stored = P_can (1 + m5 z^5 + ..)`, `m5 = -43397/807975000`) six genus-4 edge graphs are singular (graphs 1, 2 at
`[u^4]`, graphs 3, 4, 7, 8 at `[u^5]`) while the total's polar part is unchanged; per unit `c_5`: vertex `-1/8400000`, bridge (1,3)
`-1/1008000`, bridge (2,2) `1/900000`, loop 0, sum 0 (MEC s.2, s.6; S01 E12, C1, C11). Status VERIFIED (two codes, MEC s.5). Consequence
drawn by the sources: no argument about total potentials can prove Conjecture N (MEC s.2; S01 s.5).
