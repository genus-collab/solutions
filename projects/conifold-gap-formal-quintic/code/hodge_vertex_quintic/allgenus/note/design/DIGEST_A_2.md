# DIGEST_A, part 2: s.1 Definitions D11-D20

Abbreviations and labels: `DIGEST_A.md` header.

**D11. Conifold frames and the frame map.**
* PF4 in `u` (doc:fq_conifold): `p_4 D^4 + p_3 D^3 + p_2 D^2 + p_1 D + p_0` with `p_4 = u`, `p_3 = -2(1-u)`, `p_2 = -(7/5)(1-u)`,
  `p_1 = -(2/5)(1-u)`, `p_0 = -(24/625)(1-u)`. Local exponents at `u = 0`: `{0, 1, 1, 2}`; indicial polynomial `n (n-1)^2 (n-2)`
  ("checked", doc:fq_conifold; FQ2 s.1).
* Analytic solutions `f_k = u^k + O(u^3)`, `k = 0, 1, 2`: `f_0 = 1 + (2/625) u^3 + ..`, `f_1 = u - (19/60) u^3 + ..`,
  `f_2 = u^2 + (37/30) u^3 + ..` (FQ2 s.1). Logarithmic solution `g log u + h`, `h` analytic with zero `u^0, u^1, u^2` coefficients.
* Vanishing period `g = f_1 + gamma f_2`, `gamma = 7/10` "forced by the logarithm's u^-1 obstruction" (FQ2 s.1; doc:fq_conifold);
  `g = u + (7/10) u^2 + (41/75) u^3 + (1133/2500) u^4 + ..`. MS5 s.3.4: `gamma = 7/10` is now derived by hand.
* **Frame** `(alpha, beta)`: `varpi_0^cf := f_0 + alpha g + beta f_2`; `t_c := g / varpi_0^cf`. The Lagrange concomitant `B` does not
  fix `(alpha, beta)` or the scale of `g`, because `B(f_a, g) = 0` for every analytic `f_a` (doc:fq_conifold "pairing_vacuous";
  FQ2 s.1; proved in MS5 s.3.4 by the automorphism `log u -> log u + 1`). Defaults `alpha = beta = 0`, `g = u + O(u^2)`.
  Frame `(0, 0)`: `t_c = u + (7/10) u^2 + (41/75) u^3 + (9/20) u^4 + ..` (FQ2 s.1). Frame `(1/3, -2/7)`:
  `t_c = u + (11/30) u^2 + (751/1575) u^3 + ..` (NF s.6, second frame). `w_1 := [u^1] varpi_0^cf` (NF s.1 item 4; equals `alpha`);
  `gamma_t := [u^2] t_c` with `gamma_t + w_1 = 7/10` in every frame (S01 s.2 "Gauge").
* Scope widened by (CYC): any power-series solution `varpi_0` not proportional to `g` gives the cyclic factorisation of PF4;
  `varpi_0(0) != 0` is used only downstream, for the unit shape that Theorem U needs (MS5 s.3.4). Some checks used the label
  `y_0 = f_0 + a f_1 + b f_2`, not the programme's `f_0 + alpha g + beta f_2` (MS5 s.3.7; harmless).
* **Frame map Cf** (doc:fq_conifold; FQ2 s.1): the differential-ring map from N to Laurent series in `s`,
  `c_j -> D^j varpi_0^cf`, `i_j -> D^j g - j D^(j-1) varpi_0^cf - Lam D^j varpi_0^cf` (from `I_1 = C_0 Lam + Itilde_1`), `L -> 1/s`,
  `u = s^5`, `Lam ->` any value (Lam-free elements have Lam-independent images); `(C_0, I_1) -> (varpi_0^cf, g)`. Cf preserves the
  Wronskian relation (FQ2 s.1).
* **Regular at the conifold point** (FQ2 s.1): an element of N, e.g. a coloured graph term, is regular if its Cf image has no negative
  power of `s`. In the later pipeline the rows are built directly in the frame, `C0 = varpi_0`, `C1 = D t_c` (doc:ref_rows).

**D12. The hatted theory** (OU2 s.1; OU3 s.1; NF header).
* `c = (0, -1, 3, 2, 1)`; `Hhat^a := s^(1 - c_a) H^a`; `z = s^4 zhat`; `Qhat[a][i][n] := Q[a][i][n] s^(4n + 1 - c_a)`. `eta` is preserved
  on the eta-pairs `(0,3), (1,2), (4,4)`, since `c_a + c_b = 2` there.
* TFT: `Deltahat_i^(1/2) = zeta^(4i) xhat`, `xhat := x/s = varpi_0^cf(u)`. Unit `1hat = s 1` (components `Deltahat_k^(-1/2) =
  s Delta_k^(-1/2)`): the CohFT rescaling `Omega_{g,n} -> s^-(2g-2+n) Omega_{g,n}` (OU2 s.1). OU3 s.1 writes `xhat = 1 + O(u^3)`,
  which holds in the frame `(0, 0)` only (digest remark: `varpi_0^cf = 1 + alpha u + ..`).
* **(H), purity**: the hatted rows are power series in `u` (only `s^(5k)` coefficients), of valuation 0 with non-zero constant term
  (OU2 s.1: VERIFIED for all 175 rows `n <= 6`; NFS item 7; NF s.3). Final status: a consequence of `D_0 = 1` (D20).
* **Conifold symbol** `R^ := Qhat|_{u=0}`: constant, symplectic, Z/5-covariant, phase `zeta^(i(a+1))` in the normalised basis (FQ2 s.9).
* `Psi(u)`: "the constant orthonormal change of basis flat -> colour", `Psi eta^-1 Psi^T = eta^-1` (OU3 s.3); "the hatted rows are
  `Psi(u) R~(zhat; u)` with `Psi(u)` orthogonal for every `u` (it drops out of every edge)" (OU3 s.6d). [Digest reading of doc:ref_rows:
  `Psi_{aj}(u) = Qhat_{a,0}^(j)(u)`, the `zhat^0` coefficient of the rows; to be confirmed, s.7 item 7.]
* **Teleman (colour-colour) form**: `R~_{ji}(zhat; u) = (1/5) Sum_{a+b = 3 mod 5} Qhat_{b,0}^(j) Qhat_a^(i)(zhat)` (doc:ref_rows), so
  `Qhat = Psi R~`. **`M_0 := R~(zhat; 0)`**, the conifold symbol in Teleman form (NF header).
* Properties of `M_0` (FQ2 s.9): not diagonal in the colours (Z/5-circulant with phase `zeta^-n`, non-zero antisymmetric even part);
  identity row `f = 1 - (3/20) zhat + (9/800) zhat^2 + (11571/80000) zhat^3 + ..`; neither the diagonal entries nor `f` are `exp(odd)`;
  not a QRR/Gamma constant of any weight multiset (agrees with `h_{wt_0}` through `zhat^2`, departs by `37/250` at `zhat^3`).
* Clash: doc:fq_rows.WKB uses `R~` for the matrix normalised by `psi(q=0) = 1` ("R = R~ h_{wt_i}"), not for the Teleman form (s.4 C5).

**D13. Factorisations, the vertex CohFT, the outer graph sum.**
* For a symplectic `M(zhat; u)`: `R = R' M`, vertex CohFT `M.omega^` (`omega^` the hatted TFT; colour-mixing unless `M` is diagonal),
  outer graph sum `R'.(M.omega^)` over the leg-free stable graphs of genus g (NF s.5; OU3 s.2). Composition law
  `R'.(M.omega^) = (R' M).omega^`; OU3 s.2 validates it at genus 3 "for random full order-6 `R' = M_1 M_2 M_2^-1`, i.e. the outer sum
  with rows `M_1 M_2` and vertices `M_2.omega^` equals the inner tensor" of `(M_1 M_2).omega^` at `(g, n) = (3, 0)`; quoted in the proof
  as arXiv:1303.1043 Prop. 2.4 (ST2 s.2).
* Vertex classes: legs `M^-1_k`, inner edges `Ehat[(a, b)]`, leaves `That[e]` (OU3 s.2, the genus-3 vertex code).
* Edge splitting (NF s.5): `E_R = E_M + E_out`, `E_out = M^-1 E_R' M^-T = E_R - E_M`; for `M = M_0(tau zhat)`,
  `E_M(a, b) = tau^(a+b+1) Eh_0(a, b)`; the leaves of `R` are unchanged.
* Outer factors for `R' = Psi G` (UR s.6.1): edges `(1 - G(-psi)^T G(-psi'))/(psi + psi')`, leaves `z (1 - G^-1 1)` (`Psi` drops out of
  the edges by orthogonality, OU3 s.3).
* **Outer graph**: an uncoloured stable graph of the outer sum; its term contracts the full vertex tensors and so contains the sum
  over all colourings (FQ2 s.8: "the outer sum over all colourings of each stable graph"). **Vertex graph**: the single-vertex genus-g
  graph together with its R'-leaf decorations, "not the leafless class alone" (FQ2 s.10; at genus 2 the leafless class has `1/t_c`
  coefficient `+37/60000`).

**D14. The natural variable and the model point.** `tau(u) := 1 / 2F1(1, 4/5; 9/5; u)`, `1/tau = Sum_{n >= 0} 4 u^n / (5n + 4)`
(NF s.0); `tau(0) = 1`, `tau_1 = -4/9` (NF s.1 item 2); `tau(u) zhat = -(5/4) zeta^i z / (u_i - u_i^c)` (NF s.0; D5).
`v(u) := u tau^(-5/4)` (NF s.0, s.1 item 4); unique given `v(0) = 0` only up to a fourth root of unity (NFS s.4). NFS s.2 item 7: the
orientation of the integral in one formula for `tau` is opposite to the notes' sign; immaterial.

**D15. The natural form.** `M_nat := M_0(tau(u) zhat)`; `R' := R M_nat^-1 = Psi(u) G(u, zhat)`, `G = 1 + Sum_{n >= 1} u^n G_n(zhat)`,
`G_n` a polynomial in `zhat` (NF s.0, s.1 item 1). In UR s.4 (Corollary F, search only): `G = Psi(u)^-1 Ghat(u, zhat) Psi(0)`, with
`Ghat` the hatted gauge (UR s.0 line 37: `Ghat_ab = tau^((1-c_b)/4) Sum_k zhat^k u^(mu(b,k)) p_{k,b}(u) delta_{a,b+k}`).
Vertex CohFT `V_tau = M_0(tau zhat).omega^`: "the model theory in its natural variable (degree-d classes of the frozen theory
`M_0.omega^` times `tau^d`)" (UR s.6.1). Frozen theory `M_0.omega^`: TFT `Deltahat_i^(1/2) = zeta^(4i)` (FQ2 s.9). `R'` is "the hatted
form of a gauge that is regular at the conifold point in the flat frame and is computable from the periods alone" (NF s.0).

**D16. The Landau-Ginzburg model** (NF s.1 item 3; NFS items 4-6; doc:ms_rows).
* Superpotential `lambda = x^4 + u_LG/x`, primitive form `dx`; `J_k` = the normalised formal stationary-phase series of
  `Int exp(lambda/z) x^k dx`.
* Symbol: `R^_a^(i)(zhat) = N_a zeta^(i(a+1)) J_{c_a}(-4 zeta^-i zhat)`, `N = (1, -1, -1, 1, 1)`, `c = (0, -1, 3, 2, 1)`.
* Dictionary `1, H, H^2, H^3, H^4 <-> dx, -dx/x, -x^3 dx, x^2 dx, x dx` (`k = c_a`, in the window `-1..3`; NFS s.2 item 1: `x^-a dx mod 5`
  is not an identity of classes, only `x^3 dx = v x^-2 dx` holds).
* Gauss-Manin system `z d/dv Y_b = beta_b Y_{b+1}` (`b` in Z/5), `beta = (1, -1/v, 1, -1, -1)`; scalar equation
  `(v (z d/dv)^5 + 3 z (z d/dv)^4 + 1) Y_0 = 0`, local exponents `0, 1, 1, 2, 3`.
* Closed recursion for `g^a = J_{c_a}(-4 zhat)`: `g^(a+1)_n - g^a_n = ((4n - 3 - c_a)/5) g^a_(n-1)`, `Sum_a (4n + 1 - c_a) g^a_n = 0`
  (stored over Q to `zhat^40`, doc:ms_rows).
* `dx/x` (the mirror of `P(1,4)`) has the same `M_0` but another unit (NF s.1 item 3).
* "Not of Hodge/Gamma type; Mellin-Barnes kernel `Gamma(s) Gamma((s+1)/4)` in place of `Gamma(s)`" (NF s.1 item 3): context only
  ("quoted, inferred, analogy", NFS s.4).
* Frobenius manifold: the Hurwitz space of `p^4 + a_2 p^2 + a_3 p + a_4 + a_5/(p - a_1)` with form `dp`, on its residue line; flat
  metric, unit `d/da_4`, `c_555 = 1/a_5`; `H` is the residue direction; charge `d = 1/2` (NFS item 6: VERIFIED locally; NF s.0).
  Identification with `H_{0;1,4}` of arXiv:2211.12259 and its KP theorem: QUOTED, not used (NFS item 6; NF s.2).
* Normalisations in use: `x^4 + v/x` (NF), `-(x^4/4 + v/x)` (UR s.6.5), `W = p^4/4 + v/p` (S01 s.6), `x = p^4 + 4/p` (S01 C2). Only
  `eps = -4 zeta^-i zhat` is fixed by the data; the rest is the scaling freedom of the model (NFS s.2 item 7). See s.4 C9.

**D17. The model spectral curve.** `x = W = p^4/4 + v/p`, `y = p` (S01 s.6 "Route"; ST2 s.2 LD box). Topological recursion is used
only on this curve (binding rule 8; MS2 s.5 item 1: "TR on the quintic's own mirror curve stays excluded").

**D18. Regularity of an outer graph term** (OU2 s.1, OU3 s.1, NF s.5-s.6, ST2 s.2 step 6). For a leg-free genus-g outer graph
`Gamma`, the genuine term `T_Gamma` and the hatted output `G_Gamma(u) := s^(10(g-1)) T_Gamma = u^(2g-2) T_Gamma` ("in the engine's
normalisation, `G_Gamma = u^(2g-2) F_Gamma`", ST2 s.2 step 6). `T_Gamma` is **regular at `u = 0`** iff `T_Gamma` lies in `Q[[u]]`,
i.e. iff `[u^k] G_Gamma = 0` for `k = 0, .., 2g-3` (OU2 s.1 at genus 2; NF s.5: "a term is `u^-6 G(u)`" at genus 4; NF s.6 `u^-8 G(u)`
at genus 5). Polar part in `t_c` by series reversion, e.g. genus 2: `T = a_0 t_c^-2 + (a_1 + 2 c_1 a_0) t_c^-1 + ..` with
`t_c = u (1 + c_1 u + c_2 u^2 + ..)`, `a_k = [u^k] G` (OU2 s.1; genus 3 formulas OU3 s.1).

**D19. Objects of the proof chain that the statement names** (definitions belong to other digests). Cyclic systems
`z d/du Y_a = b_a Y_(a+1)` "with unit coefficients except for one simple pole" (NF s.1 item 4); Theorem U's gauge `P(u, z)`, regular
at `u = 0`, unique up to `c(z^5)` (NFS item 2); canonical gauge `P_can := P P_11(0, z)^-1` (ST2 s.2 step 2; MEC: "the Theorem-U gauge
with `P_11(0, z) = 1`"); `Y^C`, `Y^C_Q` the pure-`1/z` (vanishing-cycle) solutions of the model and of the quintic; `Pi_i` the shifted
vanishing covectors; `K = {Pi_i = 0 for i >= 1}`, `K_0 = K cap ker Pi_0` (MEC "Notation"; S01 C9 renames the slice `K_vc`, since UR uses
`K = (1, uH, H^2, H^3, H^4)` for a frame); `Omega^mod_v` the model CohFT at `v` (ST2 s.2 step 5); ML, the Model Lemma (MEC).

**D20. Normalisation objects.**
* `Q D = P_can Q^m(v(u))` with `D` a constant diagonal series, `D_i(z) = D_0(zeta^-i z)` (MS4 s.1; ST2 s.2, citing UR s.4).
* **(F_g)**: `D_0 = 1 mod z^(3g-2)`; equivalently (H) and `c^act = 1` through `z^(3g-3)`; at genus g only `P` through `z^(3g-3)` enters
  (MS2 s.1; ST2 s.2). **(N)**: `c^act = 1`, with `c^act = exp(a_1 z^5 + a_3 z^15 + ..)` the Theorem-U freedom of the actual gauge
  (S01 E4).
* Connection form (MS4 s.3): `D = 1` iff `ph^can_n(1) = h_n` for all `n`; for `5 not| n` this is (H), for `n = 5m` it is (N), with
  `f_m = [u^(4m)] ph_(5m)` (the value at the orbifold point) fixed at the conifold by `K_m = 0`.
* **(Gamma)**: `c(z; L = 1) = h(z)`, where `c(z; L)` is the full-Stirling saddle series of the Mellin-Barnes summand of the I-function
  (MS4 s.1; MS5 s.1). Final status: `D_0 = 1` to all orders, PROVED modulo (U) and (QRR) (MS5 s.1-s.2; other digests).
