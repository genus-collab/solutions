# DIGEST_A, part 5: s.4 Conventions and s.5 Evidence

Abbreviations: `DIGEST_A.md` header.

## s.4 Conventions, sign dictionaries, normalisations (as the sources fix them; differences flagged)

**C1. Weights.** Tangent weights at `p_i`: `lambda_i - lambda_j`; twist weight `-5 lambda_i` ("numerator twist enters with a MINUS
sign", doc:fq_hodge.wt; FQ2 s.1). General rule: "an Euler-class twist contributes minus its fibre weight, an inverse Euler class plus
its fibre weight" (MS5 s.3.5). The sign is pinned by `h_1 = ph_1 = -3/20` (MS5 s.2; MS4 s.2: "pinned by (H) at n = 1"). No conflict.

**C2. Pairing.** `eta(H^a, H^b) = 5 delta_{a+b = 3 mod 5}`, `eta^{ab} = (1/5) delta` (FQ2 s.1; doc:fq_rows). In the basis
`h^a = N_a H^a` of doc:ms_rows the same (`N_a N_b = 1` on every eta-pair). In the hatted basis `eta` is preserved because
`c_a + c_b = 2` on the eta-pairs (OU2 s.1).

**C3. Square roots.** `zeta^(1/2) := zeta^3`, so `Delta_i^(1/2) = zeta^(4i) C_0/L` and `(W_i^0)^(-1/2) = zeta^i` (doc:fq_rows.Rows).
This fixes the phases `zeta^(ia)` of `P` and `zeta^(i(a+1))` of the normalised rows and of the symbol (FQ2 s.9; D16).

**C4. Row normalisations (flag).** Five objects: `R_{ai}` (flat component `a`, colour `i`, f-basis), `P = Delta^(1/2) R`, `X_a`
(`P = (C_0/L) X`), `Q = P/Delta^(1/2)` (doc:ref_rows), `Qhat = Q s^(4n+1-c_a)` (OU2 s.1). Valuations: `val_s P_{a,n} = -4n + c_a`
(FQ2 s.3), `val_s Q_{a,n} = -4n + c_a - 1` (OU2 s.1), consistent with `val_s Delta^(1/2) = 1`. OU2 s.1 also writes
`Q[a][i][n] = [z^n] eta(H^a, R(z) e~_i)`; read literally this lowers the index `a` with `eta`, which does not match doc:fq_rows'
`R_{ai}` without a relabelling `a <-> 3 - a` and a factor 5. Not resolvable from the sources read: s.7 item 6.

**C5. The letter `R~`.** Teleman colour-colour form `R~(zhat; u)` (NF, doc:ref_rows) versus the I-function-normalised matrix with
`psi(q=0) = 1` in doc:fq_rows.WKB ("R = R~ h_{wt_i}"). The note needs two symbols.

**C6. Edges.** Point frame, colour basis: `[delta_ij Delta_i - Sum P_a^(i)(-z) eta^{ab} P_b^(j)(-w)]/(z + w)` (doc:fq_rows.edge_coeffs,
citing a draft formula (1.1)); model, flat basis: `(eta^-1 - R^-1(psi) eta^-1 R^-1(psi')^T)/(psi + psi')` (doc:ms_rows); outer:
`(1 - G(-psi)^T G(-psi'))/(psi + psi')` (UR s.6.1, `Psi` dropped). The quoted convention is arXiv:1303.1043 Def. 2.2 (ST2 s.2). Digest
remark: the three agree through symplecticity; the writers should state the PPZ form once and derive the colour-basis form.

**C7. Leaves.** Point frame `Dil_i(d) = (-1)^d [z^(d-1)] P_0^(i)`, `d >= 2` (doc:fq_rows.Dil); model `T(e) = -R^-1_(e-1) h^0 =
(-1)^e g^0_(e-1) h^(1-e)`, `e >= 2` (doc:ms_rows); outer `z(1 - G^-1 1)` (UR s.6.1). All are the translation `z(1 - R^-1(z) 1)` in different
bases (digest remark; quoted convention arXiv:1303.1043 Def. 2.5, ST2 s.2).

**C8. Hatted grading.** `c = (0, -1, 3, 2, 1)` ("the hatted weights are the Euler grading of the model", NF s.1 item 3);
`e = c - 1 = (-1, -2, 2, 1, 0)` (doc:ref_rows); `1hat = s 1`; `Omega_{g,n} -> s^-(2g-2+n) Omega_{g,n}`; `xhat = varpi_0^cf(u)` (OU2 s.1).

**C9. Model normalisation (flag).** Four forms (D16). "u_LG = 4u, z_LG = -4z exactly" (LG) and "c^4 = -1/4" (UR and referee) "are the
same statement; only `eps = -4 zeta^-i zhat` is fixed by the data, the rest is the scaling freedom of the model" (NFS s.2 item 7). The TR
route uses `x = W = p^4/4 + v/p`, `y = p`, `d/dx = p^2/(p^5 - v) d/dp`, legs `Phi_0^b = p, log p, v/p, p^3, p^2` (ST2 s.2 LD box).
Proposal: one normalisation in NOTATION.md (the TR one), with the dictionary to `J_k` stated once.

**C10. `H` and `e` (adjudicated).** The project's `H` is `-d/dv` at its frozen point `t = -1 <-> v = -1`; hence `<H>_1 = -1/12`,
`<e>_1 = +1/8` (tr, model, xyswap, theory) against `+1/12`, `-1/8` (UR s.6.3 and the plan). "The difference is convention"; advice: fix
`H = d/dv` and `e = tau_1(1) + v tau_0(H)` (S01 C3; MEC s.6; ST2 header adopts it).

**C11. Euler-dilaton block sign.** Flat `-pi z (uH + z) (x) (uH - z)` versus hatted `alpha zhat (H - zhat) (x) (H + zhat)`: the
difference is that "the Fourier vector of H is minus the flat class (`N_1 = -1`)", not an ordering convention (NFS s.2 item 8).

**C12. Sign of `c_5` (adjudicated).** One worker note used both signs; unified as `P_stored = P_can (1 + m5 z^5 + ..)`,
`m5 = -43397/807975000` (S01 C11; MEC s.6). NFS s.5 prints `c_5 = 43397/807975000` and MEC s.2 "c_5 = -43397/807975000 away from
canonical": both refer to the same number under the two conventions; use S01 C11.

**C13. Vertex constant.** Four equal forms (S01 C2; A16). The `5^(1-g)` is the TR-to-frozen factor `(-5)^(1-h)` with sign, for `x = W`
(S01 C2; ST2 s.2 step 7.1: `<q_1..q_n>_h = (-5)^(1-h) <Phi_1..Phi_n>^TR_h`).

**C14. TR conventions (for the TR digests).** The kernel is minus the kernel printed in arXiv:math-ph/0702045 (4-3); `F_2^TR = -1/240`
on `x = W` (ST2 s.1 table M8); string-equation sign `+d/dx` (S01 C8); `x = -W` flips `psi -> -psi` only up to `(-1)^(3g-3+n)`, and the
metric sign is `(-1)^(g-1+D)` (S01 s.2 "tr"). LD constants `c = (-1, 1, -1, -3, -2)` with `q_{Phi_k^b} = -c_b tau_k(H^b)` (MS2 s.2).

**C15. Other adjudicated conventions (other digests).** `lambda` in Cor. U: define `r = -1/lambda^2`, `lambda = f(0)|g(0)|` (S01 C10).
Vanishing covector: `PiC(tau_k H^b) = (-1)^(k+1) [z^-k] Y^C_b`, closed form `1/(4^q q! (m-q)!)`, support `b = 1 - m` (S01 C4).

**C16. Frame labels.** `(alpha, beta)` always means `varpi_0^cf = f_0 + alpha g + beta f_2`; checks labelled `y_0 = f_0 + a f_1 + b f_2`
(MS5 s.3.7) and the block coefficient called `alpha` in NF s.1 item 4 and NFS item 11 are different objects.

**C17. Scale of the vanishing period.** `g = u + O(u^2)`; the scale is free (doc:fq_conifold), and `c_g` is stated for this scale (FQ2 s.4).

**C18. Hatted versus flat.** `G_Gamma = u^(2g-2) F_Gamma` (ST2 s.2 step 6). Theorem S15 uses "Lemma 3's description of the hatted
variables", which "matches the fq_sser/fq3_msum docstrings" (ST2 s.2 step 6): the note's D12 must be that description.

**C19. Rationality.** Per coloured term the coefficients lie in `Q(zeta)`; only Z/5-orbit sums are rational (FQ2 s.3); `Q(zeta)` entries
are recovered from Galois conjugates `sigma_m` = colour relabelling `i -> m i` (OU2 header). In common-weight frames the per-colour
series are not u-series, so per-colour polar parts are defined only in the covariant frames (FQ2 s.3). The theorem is per outer graph.

**C20. Status vocabulary.** NF uses ARGUED for "written and audited"; NFS uses PROVED-mod(X); MEC onward use PROVED. Final statuses are
those of MS5; NFS s.4 and S01 s.2 ("Record notes now outdated") list the upgrades of NF s.3 (model gap for all g; bookkeeping derived;
normalisation verified through `z^24`).

**C21. Letter clashes in the sources** (rename in NOTATION.md).

| letter | meanings found |
|---|---|
| `s` | `u^(1/5)`; `dT_P - v` (S01 E5); the series `s(z)` (MS4 s.3); Hodge coefficients `s_k` (doc:fq_hodge) |
| `tau` | natural variable `tau(u)`; descendants `tau_k(H^b)`; mirror-map value `tau_0` (FQ2 s.1); KP tau function |
| `alpha` | frame parameter; Euler-dilaton block coefficient (NF s.1 item 4); `alpha0` in LMI (MS5 s.3.1); exponents in LD |
| `c` | hatted weights `c_a`; LD constants `c_b`; `c_g`; Theorem-U freedom `c(z)`; saddle series `c(z; L)`, `c_n(L)`; `c^act`; `c_j = D^j C_0` |
| `K` | the frame `(1, uH, H^2, H^3, H^4)` (UR); the slice `K` (MEC; renamed `K_vc`, S01 C9); `K_q` (MS4 s.3); window size `K = 8` (NF s.6) |
| `G` | gauge factor `G(u, zhat)`; series `G_Gamma(u)`; generating function `G` (S01 s.6); `G = log(..)` (MS4 s.4.3); `Gh` (doc:ms_rows) |
| `M` | `M = z D` (PF); factor `M`; `M_0` |
| `L` | `(1 - 5^5 q)^(-1/5)`; leg class `L` (S01 E9); lemmas (L1)-(L3) |
| `D` | `q d/dq`; series `D(z)`; D(ii), D(iii); `d/dx` (ST2 LD box); psi-deficits `D_v` |
| `H` | hyperplane class; Hessians `H_7`, `H_vv`; hypothesis (H) |
| `P` | gauge `P(u, z)`; rows `P_a^(i)`; `P^4`; Lemma P, Theorem P(b) |
| `Q` | Novikov variable; rows `Q`; the field |
| `T` | mirror map; torus; graph term `T_Gamma`; leaves `T(e)`; `T_P`; Lemma T |
| `N` | ring N; signs `N_a`; `N_{g,d}`; Conjecture N; input (N); `N_k` |
| `W` | `W_i^0`; Wronskian denominator `W`; model superpotential `W`; "rule W" |
| names | "Theorem G" (gauge identity vs gap theorem: use (G) and Gap Theorem); "Lemma K"; "ML1/ML2" swapped; "Theorem R" also the conifold theorem (S01 C9) |

## s.5 Evidence for the appendix (finite checks; ranges, fields, counts; each with its source)

Fields: `F_p` with `p = 2^89 - 1`, `2^61 - 1` and rational reconstruction; `2^521 - 1` in one referee check; exact `Q` or `Q(zeta)` where
stated (FQ2 header; OU2 header; NF header).

**V1. Calibration** (FQ2 s.2). Graph sum against direct localisation: `F_2(0) = N_{2,0} = -19/9600`, `N_{1,1} = -305/6`,
`N_{1,2} = -535775/12`, `N_{2,1} = 2129/5760`, `N_{2,2} = 203401/144` (both sides); predictions `N_{1,3} = -1540037125/36`,
`N_{1,4} = -1082838731875/24`, `N_{2,3} = 1016639975/384`, `N_{2,4} = 181712650375/48`; `[Q^0] d_T F_1 = -3/8`; `F_3(0) = 1/224000`
(localisation); genus-0 numbers independent of the weights at four choices; the `d = 1` genus split per graph (vertex graph `-24161/5760`,
two genus-1 vertices `2629/576`); controls that must fail fail. `N_(2,0)`, `N_(3,0)` also from arXiv:1803.01409 eq. (26) (MS5 s.4.1).
Not done: `F_3(0)` from the project's R-matrix graph sum (MS5 s.5 item 3).

**V2. Genus 2, point frame, per graph** (FQ2 s.3; OU2 s.1; graphs 1 vertex, 2 genus-1 loop, 3 genus-0 two loops, 4 g1-g1, 5 g0 loop -
g1, 6 theta, 7 dumbbell). `t_c^-2`: `1807/3072000, -811/256000, -57/64000, 843/1024000, 3/51200, 87/32000, 9/12800` (sum `1/1200`);
`t_c^-1`: `-791/15360000, -1771/3840000, 41/12800, 743/15360000, -41/153600, 37/160000, -173/64000` (sum 0). 7 stable graphs, 115
coloured graphs; identical tables at both primes; mixed colourings `Q(zeta)`-valued (graph 4 at (0,1): `1/40000 - zeta^2/30000 -
zeta^3/30000`).

**V3. Genus-2 total** (FQ2 s.4). `F_2^cf = (1/1200) s^-10 - (7/6000) s^-5 + 97/72000 - (6809/18000000) s^5 - .. = (1/1200) t_c^-2 +
0 t_c^-1 + 31/30000 + O(t_c)`; `[s^-5]/[s^-10] = -7/5 = -2 gamma`; frames `(0,0), (1,0), (0,1), (2/3,-1/7)`: `c = 1/1200`, `c' = 0`,
constant terms `31/30000, 31/30000, 2749/120000, -251/120000`.

**V4. Negative results at genus 2** (FQ2 s.3, s.6-s.11): A18; 110 coloured edge graphs; constant covariant scalar, diagonal (10 constants)
and symplectic (25 parameters) classes; both primes. Covariant-scalar O(u) solution `(b_1, b_3) = (113/75, -21/125)` (left placement
`(-113/75, 2/15)`), inconsistent with genus 3 (FQ2 s.11).

**V5. Natural form, genus 2 and 3** (NF s.1 item 5; UR s.6.2; NFS s.0). `M` mod `(zhat^7, u^4)`; frame `(0,0)` at both primes, frame
`(1/3, -2/7)` at `2^89 - 1`. Genus 3 (42 graphs; 1, 2, 5, 9, 12, 8, 5 with `E = 0..6`, OU3 s.2): all 41 edge graphs `O(u^4)`; vertex
`u^-4 (1/25200 - u/9000 + 407u^2/3780000 - 151u^3/3780000) = (1/25200) t_c^-4 + 0 + 0 + 0`; second frame
`u^-4 (1/25200 - 11u/189000 - ..)`, same polar part. Genus 2: six edge graphs `O(u^2)`, `[u^2] = 479/756000, 1/2250, -127/2016000, 0, 0, 0`;
vertex `[u^0] = 1/1200`. Non-kinematic content: 9 rational numbers per frame (NFS item 14).
