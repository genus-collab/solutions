# DIGEST_A, part 6: s.5 Evidence (continued: V6-V13) and s.6 Literature

Abbreviations: `DIGEST_A.md` header. V1-V5: `DIGEST_A_5.md`.

## s.5 (continued)

**V6. Genus 4** (NF s.5). 379 leg-free graphs (1, 3, 7, 21, 43, 75, 89, 81, 42, 17 with `E = 0..9`); a term is `u^-6 G(u)`. Method:
edge splitting `E_R = E_M + E_out` (D13): the point-frame sum over graphs with every edge marked inner or outer, regrouped by contracting
the inner edges, is the outer sum `R'.(M.omega^)`; all vertices are Witten-Kontsevich numbers with one colour; `M` exact in `u`, mod
`zhat^(3g-2)`; series mod `u^(2g-2)`. Frames `(0,0)` and `(1/3, -2/7)`, both primes. Total, frame `(0,0)`:
`u^-6 (1/180000 - 7u/300000 + 701u^2/18000000 - 1459u^3/45000000 + 923639u^4/67500000000 - 171583u^5/67500000000)`
`= (1/180000) t_c^-6 + 0 + .. + 0`; second frame `u^-6 (1/180000 - 11u/900000 - 79u^2/378000000 + ..)`, same polar part. All 378 edge
graphs `O(u^6)` in both frames at both primes; the vertex graph carries the polar part; `1/180000 = (-B_8/48)/125`.
Content: 229 graphs with `E >= 6` regular a priori (A9); of the 149 with `E <= 5`, 81 (3, 19, 59 with 3, 4, 5 edges) receive only pieces
that vanish identically (a selection rule of the true `E_out`; a generic `E_out` of the same valuation gives non-zero pieces) and 68
(3, 7, 18, 24, 16 with 1..5 edges) vanish by cancellation among 6767 non-zero pieces. Controls (`2^89 - 1`, frame `(0,0)`): bare symbol
`M = M_0` (`tau = 1`): all 149 graphs with `E <= 5` singular; natural `tau` with only `tau_5` shifted by 1: exactly the three one-edge
graphs fail, at `[u^5]`, and the vertex graph acquires a `t_c^-1` term. Independent checks (NF s.5): a second code sharing no module
reproduces, at both primes in frame `(0,0)`, all 379 outer values, all 379 point-frame terms, every non-zero piece and both controls
(not the second frame); 26 code mutants, all caught (23 by the genus-2/3 regression, `tau_4`, `tau_5` and powers of `tau` at `zhat >= 7`
by the genus-4 checks). Caveat (NF s.5): the identity "regrouped split sum = outer graph sum" is an argument, checked against genuine
vertex tensors only at genus 2 and 3.

**V7. Genus 5** (NF s.6). 4555 graphs (1, 3, 11, 34, 100, 239, 492, 784, 1002, 926, 632, 260, 71 for `E = 0..12`); a term is
`u^-8 G(u)`. Frame `(0,0)`, both primes; frame `(1/3, -2/7)`: regularity at `2^89 - 1` only, total at two primes and certified.
Total (frame `(0,0)`, by CRT from both primes): `u^-8 (1/660000 - 7u/825000 + 199u^2/9900000 - 717u^3/27500000 + 449357u^4/22500000000 -
5593217u^5/618750000000 + 14686042363u^6/6496875000000000 - 8370858511u^7/32484375000000000) = (1/660000) t_c^-8 + 0 + .. + 0`;
`1/660000 = |B_10|/(10 . 8 . 5^4)`. Second frame total: `u^-8 (1/660000 - u/225000 + 323u^2/207900000 + 66841u^3/15592500000 - ..)`, which
differs at every order `u^1..u^7` from the first while the polar part is again `(1/660000) t_c^-8`. All 4554 edge graphs regular.
Content: 2891 graphs with `E >= 8` regular a priori (A9; every entry of `E_out` has valuation exactly 1, 1950 entries); 1036 (7, 71, 288,
670 at `E = 4..7`) receive no non-zero piece (selection rule); 627 (3, 11, 34, 93, 168, 204, 114 with `E = 1..7`) vanish by cancellation;
3142 polar coefficients at stake, 1655 of them in rows that receive pieces: 517 receive no contribution at that order, 1138 are genuine
cancellations (pieces per cancelling row: min 31, median 347, max 4525; the one-edge rows are sums of 3479, 2296, 4525 pieces). The
selection rule is frame independent on rows (same 628 rows, same 627 zeros in both frames) and frame dependent on pieces (190 pieces
differ, all with `E(t) >= 7`). Controls (restricted to outer graphs with at most two edges): bare symbol: 14 singular = every graph
computed, vertex spoilt at `t_c^-7..t_c^-1`; `tau_7` shifted by 1: the three one-edge graphs singular at `[u^7]` only, vertex spoilt at
`t_c^-1` only. Control analysis: a shift of `tau` at order `k` moves an `E`-edge row at valuation at least `k + E - 1`, visible iff
`k + E - 1 <= 2g - 3`; `[u^1](M_0 - nat) = (4/9) [u^7](nat_t7 - nat)` on the vertex and the three one-edge rows; hence within
`M = M_0(tau zhat)`, `tau` is the unique rescaling mod `u^8` for which the one-edge graphs are regular. Independent checks (NF s.6): own
CRT and lattice reconstruction, own Frobenius solution and own `t_c`; graph list complete and non-redundant by an independent induction;
`|Aut|` recomputed; a module-disjoint second implementation on 1474 of 4555 inner graphs (93551 pieces at `2^89 - 1`, 19591 at
`2^61 - 1`, 0 mismatches); a third prime `2^521 - 1` on 31 inner graphs; per-inner-graph additivity for all 4555 inner graphs at both
primes and in the second frame (13665 identities). Caveats (NF s.6, to be stated plainly in the appendix): the arithmetic evidence is
probabilistic (heights up to `2^136` against `p89 p61 = 2^150`; failure probability about `2^-138`); the second frame's regularity rests
on one prime; `tau_6` is untested by any control at any genus, and so are the orders `zhat^10..zhat^12` of `M_0`; the two controls
discriminate on 14 of the 1663 content graphs; the regrouping identity is the genus-4 caveat of V6; the graph list is shared.

**V8. Factorisation** (NF s.1 item 1; NFS item 1). Project rows: `zhat <= 6`, `u <= 4`, both primes, two frames (11 vanishing matrix
coefficients, 0 free parameters); rebuilt at `zhat <= 16`, `u <= 11` (`2^89 - 1`) and `zhat <= 12`, `u <= 8` (`2^61 - 1`). Independent
rows: `zhat <= 16`, `u <= 11` at three primes; `zhat <= 24`, `u <= 17` at two primes (215 vanishing coefficients). Left placement false;
`tau = 1`, `(1-u)^(4/9)`, `(u/t_c)^(4/5)` fail at the first order where they differ (three primes) (NF s.1 item 2). Degrees
`floor((5n+2)/4)` (bound) VERIFIED `zhat <= 24`, `u <= 17`; seven frames at `zhat <= 12`, `u <= 8` (NFS items 10, 13).

**V9. Symbol = model** (NFS items 4-5). All rows `zhat <= 18` at two primes, identity row `zhat <= 40`, two stationary-phase engines;
only the dictionary `k = c_a` matches, both primes. Hurwitz identification exact at 4 random points, floating point at 6 (NFS item 6).

**V10. (H), the normalisation, `D_0 = 1`, (Gamma).** (H): mod `(z^17, u^12)` at `2^89 - 1` and `(z^13, u^9)` at `2^61 - 1` on the project
rows; independent rows pure by construction to `zhat^24`, `u^17`, with the QRR identity automatic for all `n <= 24`, `5 not| n` (NFS
items 7, s.2 item 4). `P_11(0,z) = 1` at `z^5, z^10` (three primes), `z^15, z^20` (two primes) (NFS item 8); the gauge computed from the
periods reproduces all 397 Fourier components of `G` for `zhat <= 14`, `u <= 11` as exact rationals with no fitted constant (NFS item 9).
`c^act = 1` in the window `zhat <= 24`, `u < 18` (`a_1 = a_3 = 0`; two primes, two frames, exact on 397 components; "covers genus <= 9",
S01 E4); through `z^35` (`a_1, a_3, a_5, a_7 = 0`) at two primes and through `z^150` (two primes, two frames) (MS4 s.4.1). `D_0 = 1` exactly
over Q through `z^204` (`g <= 69`), two codes, and through `z^125` by a third; mod two primes for `g <= 335` (MS4 s.2). (Gamma) mod
`2^127 - 1` and `2^107 - 1` through `z^500`, exactly through `z^60` (MS4 s.2, s.4.1; MS5 s.4.1); a method-independent interpolation code
through `z^30` (MS4 s.4.2); the two-variable chart against the full-Stirling series at four values of `L` and against `h` at `L = 1`,
exact over Q through `z^16` (MS5 s.4.2); `Phi^SP = Phi^can` through `w^80` at two primes in two frames; `V = C0 (Y^C_Q)_0` exactly
through `u^40` (MS5 s.3.3); (CYC) through `u^56`, `u^60` (four frames), `u^80` (MS5 s.4). Twelve theories, 40 colours: the mechanism's
predictions hold exactly to `z^24` / `z^20` and fail exactly where a hypothesis fails (MS5 s.3.5).

**V11. Per-graph sensitivity and universality.** Stored normalisation `m5 = -43397/807975000`: six genus-4 edge graphs singular, total
unchanged (A17; MEC s.2, s.5; S01 E12). Universality and B-variation: A27 (S01 E12; MS2 s.3 (c)). (G) VERIFIED over Q, `z <= 20`,
`u < 29`, two frames (S01 E1).

**V12. Model side** (pointers; other digests). D(iii) exact for weights `<= 400`, two primes `<= 1000` (MS2 s.2); log-free vectors 10400
and 16200 all in L (MS2 s.3 (a)); TR correlators `omega_(4,2)`, `omega_(3,3)`: no non-zero pairing with negative `v`-exponent among
17416 singular tuples per prime (MS2 s.3 (b)); ML on K at genus 0-4: 1395 singular K-correlators, all zero (S01 E10); the Gap Theorem
`h <= 6` (ST2 s.1 M8); flat = hatted at genus 3 on 336 coefficients per prime and frame (MS2 s.3).

**V13. Wording rules for the appendix, from the sources.** Only edge-graph regularity tests the factorisation (A26); the value `c_g`
and the gap of the total are `M`-independent (NF s.6); "regular a priori" graphs (A9) are an argument, not evidence (NF s.6).

## s.6 Literature named in the sources read

Columns: id; statement number; what is used; read in full text by the sources?; source; root refs.bib key (scan of `refGW/refs.bib`;
"none" = no entry found by id, author or title; "(by title)" = key matched by title only, eprint field absent, confirm).

| id | statement | use | full text read? | source | root key |
|---|---|---|---|---|---|
| arXiv:1702.06096 (LP) | s.0.2; eq. (5) | definition of the formal quintic; mirror map "(LP)"; the LP shape of `F_2` | eq. (5) "as read by the referees; statement numbers not re-checked" (MS4 s.2) | FQ2 s.1, s.4; MS4 s.2 | none |
| arXiv:1803.01409 (LP, formal quintic) | v2 eq. (26); s.2.4.4 eq. (18), Prop. 5 | localisation vertex class, vertex-level (QRR); `N_(2,0)`, `N_(3,0)`; I-function at `H = lambda_i` | lit worker (MS5 s.2) | MS5 s.2, s.4.1; stage4/gamma/GAMMA.md l.148 (grep) | none |
| arXiv:1509.06954 | Thm 8.1 (case (8-1), `n = 4`, `a = (0,0,0,0,v)`, `Lambda = -5 kappa`); (2-35) | (Q): Givental = TR for the model | hypotheses checked in Stages 1-2 (MS5 s.2 table) | S01 E8, s.6; ST2b s.6.2 | Dunin-Barkowski:2015caa |
| arXiv:2212.00320 | Thm 7.1 (`y = z`, rational `x`, simple critical points) | (Q): L-regularity; a gap route | as above | S01 E8-E9; ST2b s.6.2 | none |
| arXiv:math-ph/0702045 | Thm 5.1 with the linear loop equation of s.4; s.5; kernel (4-3) | (Q) on the model curve only; variational formula | as above | MS2-MS5 s.2; S01 s.6; ST2 s.1 M8 | Eynard:2007kz |
| arXiv:1303.1043 (PPZ) | Def. 2.2, Prop. 2.4, Def. 2.5, Prop. 2.9 | (Q): edge formula, composition, translations | "located in the PDF by the flat referee" (ST2 s.2) | ST2 s.2; FQ2 s.8 | none |
| arXiv:2502.18115 | Cor. 1.2, formula (1.3) (example s.3.2 has a sign slip) | optional gap route | - | ST2 s.1 M8; MS5 s.2 | none |
| arXiv:2211.12259 | Thm 2; KP theorem | ML on the `H^2` tower (not needed); identification with `H_{0;1,4}` | "read only in part, not used" (NF s.2) | NF s.2; MEC s.1; NFS item 6 | none |
| arXiv:0710.0889 (Zagier-Zinger) | Thm 4(ii); expansion (5); Thm 4(i) | degree bound of `ph_n`; `Phi^ZZ`; I-function normalisation | - | MS4 s.3; stage5/lit/LIT.md l.75 (grep) | none |
| arXiv:math/0108100 (Givental) | Thm 9.1; s.10, Lemma 10.4; v2 Thm 4.1 | Bernoulli form (untwisted; superseded as (QRR) source); chart mechanism (outline, prior art); vertex-level (QRR) | - | MS4 s.2-s.3; MS5 s.2, s.3.6 | MR1901075 (by title) |
| arXiv:1812.11908 (GJR) | v1 eq. (17), proof of Lemma 5.2; Remark 7.2; Lemma 7.7, eq. (37) | (QRR) for the formal quintic; prior art for (Gamma3) and Ginf | not in the typeset PDF (MS5 s.2 table: "human reading of GJR (17) in the typeset PDF") | MS5 s.1-s.3.6 | none |
| arXiv:1709.07392 (GJR) | v1 eqs. (19)-(20) at `t = 0` | (QRR) source | - | MS5 s.2 | none |
| arXiv:math/0008067 | Thm 1 with its Remark | cited by GJR for (QRR) | - | MS5 s.2 | none |
| arXiv:math/0110142 | Thm 1, Cor. 1 | cited by GJR for (QRR) | - | MS5 s.2 | MR2276766 (by title) |
| arXiv:0712.0160 (Teleman) | v3 Theorem (2) (s.1.9), Prop. 8.5(ii) (s.8.2) | uniqueness of R for the CohFT | "read only through a summariser" (MS5 s.2) | MS5 s.2, s.5 | MR2917177 |
| arXiv:0903.1463 (Iritani) | v3 Lemma 4.20 | prior art (exact q-expansion at `z = 1`) | - | MS5 s.3.6 | MR2553377 (by title) |
| DLMF 5.11.8 | 5.11.8 | Stirling-type expansion; a formal proof makes Ginf purely formal | - | MS5 s.3.7; task text | none |
| arXiv:2312.16950 | - | candidate route for (V); not used | - | ST2b s.6.4 | none |
| Mumford; Faber-Pandharipande | no numbers | Mumford's relation and GRR (vertex integrals, (QRR) derivation) | - | doc:fq_hodge; MS5 s.2 | MR717614; MR1728879 |
| Witten; Kontsevich | no numbers | Witten-Kontsevich numbers at the vertices | - | doc:fq_hodge | Witten:1990hr; Kontsevich:1992ti |
| Graber-Pandharipande | no numbers | localisation (calibration only) | - | FQ2 s.2 | MR1666787 |
| "HKQ" | none | conifold-frame style; "the HKQ constant" | - | FQ2 s.1, s.4 | none (s.7 item 14) |
| Chester-Friedman-Ursell | none | the CFU map in Theorem P(b) | - | MS4 s.3 | none |
| Hurwitz Frobenius manifolds | none in the sources | context for D16 | - | NF s.0; NFS item 6 | candidate Dubrovin:1994hc (hep-th/9407018), not cited by the sources |
