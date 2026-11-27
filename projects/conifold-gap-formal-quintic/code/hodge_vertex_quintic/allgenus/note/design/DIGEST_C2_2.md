# DIGEST_C2, part 2: s.2 Statements, model side (the model as TR, the lattice, log legs, ML_K)

Continues `DIGEST_C2.md` (source keys and definitions there). Proposed label prefixes follow the write-up plan of SYN2b s.6.3:
`trm` (plan s.5: the model as TR, LD, D(ii), D(iii)), `mlk` (plan s.6: L-regularity, Theorem V, Gap, ML_K), `flh` (plan s.7),
`mai` (plan s.8), `uni` (plan s.9). Each entry: label, type, statement (final corrected form), source, status in the sources,
corrections applied, dependencies. Proofs: s.3 (parts 4-5).

## s.2.1 The model as topological recursion

**`thm:trm-ident` (thm; quoted + derived) TR identification of the model CohFT (M1).** For `v != 0` and `2h - 2 + n > 0`,
`omega_{h,n} = (-5)^(h-1) Sum <prod_i tau_{d_i}(H^{a_i})>_h prod_i d xif_{a_i,d_i}`, the brackets being the ancestor correlators of the
frozen model `M_0.omega^` at `v`. Equivalently `<q_1 .. q_n>_h = (-5)^(1-h) <Phi_1 .. Phi_n>^TR_h` whenever `q_i = q_{Phi_i}` (D5).
Free energies: frozen `F_h = (-5)^(1-h) F_h^TR` (`h >= 2`), by the dilaton equation on both sides (dilaton leg = LD at `(b,k) = (0,1)`).
* Source: SYN2 s.1 M1, s.2 step 7.1; TRM s.0 item 2, s.1, s.2.2-s.2.4; LL s.1 ("Which correlators").
* Status: PROVED mod arXiv:1509.06954 Thm 8.1 (hypotheses checked, TRM s.1) + uniqueness of the homogeneous R-matrix (TRM s.2.3).
  VERIFIED (Stage 1: 72120 correlators, two primes; Stage 2: 52 frame one-point functions, `g = 2, 3`, exact over `Q(zeta5)`).
* Corrections: K1 (the phrase "the kernel of arXiv:math-ph/0702045" in TRM means the engine convention D2, not the printed (4-3));
  S01 s.7: TRM's "x = -W flips psi -> -psi" holds only up to `(-1)^(3g-3+n)`.
* Hypotheses as checked (TRM s.1): DNOPS Thm 8.1, superpotential (8-1) `f = p^n + a_2 p^(n-2) + .. + a_n + a_(n+1)/(p - a_1)`, here
  `n = 4`, `a = (0,0,0,0,v)`, `v != 0`: critical points `p^5 = v/4`, non-degenerate (`lambda'' = 20 p_c^2`), critical values `5 p_c^4`
  pairwise distinct, so semisimple with one critical point per singular fibre; the shift `p_Dub = p - a_1/(n+1)` is trivial (`a_1 = 0`).
  Metric: DNOPS (2-7) is MINUS the metric (1.6) of arXiv:2211.12259 and of LG_MODEL s.5.1; the two CohFTs differ by
  `Omega_g -> (-1)^(g-1) Omega_g`. Normalisation `Lambda = -5 kappa` (`kappa = 1` for `x = W`): factor `(-5)^(g-1)`, no psi rescaling.
* Depends on: DNOPS (2-26)-(2-35), (2-7), Thm 8.1; the `(U, V)`-recursion `[R_(m+1), U] = (m + V) R_m`; agreement of `U` (critical values)
  and of `R_1` (TRM s.2.3: obtained from a finite check, see s.7 item G2); `xif_a` of D4.

**`lem:trm-xfun` (lem) Functions of x pair to zero (M6).** For `f` a polynomial in `x` (coefficients in `Q(v)`) and `2h - 2 + n > 0`,
the pairing of `omega_{h,n}` with `f(x)` in one slot vanishes, whatever the other slots carry. Hence `L` acts modulo `Q(v)[x]`.
* Source: SYN2 M6; LL s.1. Status: PROVED: "linear loop equation (standard; arXiv:math-ph/0702045 s.4; statement number not
  re-checked)". [DW]: the printed linear loop equation is EO Theorem 4.4, (4-23)-(4-24) (p. 24). The derivation of M6 from it is not
  written in any source (s.7 item G3). Evidence: consistency check only (LR F4; SYN2b O12).

**`prop:trm-legdict` (prop) The leg dictionary LD (M2).** For all `a, b in {0..4}` and `k, d >= 0`:
`<Phi_k^b, d xif_{a,d}> = -c_b delta_ab delta_dk`, `c = (-1, 1, -1, -3, -2)`. Hence `q_{Phi_k^b} = -c_b tau_k(H^b)`; the TR leg of
`v^m tau_k(H^b)` is `-c_b^-1 v^m Phi_k^b`, and the TR leg of `s in S_w` is `-Lambda(s)`.
Consequences: `tau_0(1) <-> y = p` (the string pairing `<y, d xif_{a,d}> = delta_{(a,d),(0,0)}` of TRM s.2.2); the dilaton leg
`Int y dx = -Phi_1^0` has `q = -tau_1(1)`; the log leg `-log p` has `q = tau_0(H)` (= LL s.1 identity (X):
`<-log p, d xif_{a,d}> = delta_{(a,d),(1,0)}`); `e = tau_1(1) + v tau_0(H)` has TR leg `-p^5/5`.
* Source: SYN2 s.2 step 7.2 (boxed proof), M2; LL s.1 (X); DR s.2 Correction 1 (the `e` example).
* Status: PROVED (SYN2 s.2 step 7.2, the critic); VERIFIED (`cr2_legdict 24`, `cr2_e2e`). Not separately refereed (s.7 item R1).
* Corrections: O26 (Stage 1 assumed LD), DR Corrections 4(i), 5; K10 (LD, not XYSWAP's `tau_0(1) <-> -z` with step `+Int dx`, is the
  dictionary of record); confirms the sign of S01 E11(ii).
* Depends on: D4 (the explicit `xif_a`); integration by parts and the residue theorem on P^1.

**`lem:trm-lemG` (lem) Generating function (Lemma G).** For all `b, k`: the `log p`-coefficient of `Phi_k^b` is
`Sum_{i=0..k} R^b_(k-i) (-x)^i/i!`, with `R^b_n := Res_{p=0} W^n dPhi_0^b/n!` (a monomial in `v`). Source: DR s.3. Status: PROVED.
Depends on: `lem:trm-closure` (M is a `Q[v,x]`-module).

**`lem:trm-resid` (lem) Residue identity.** `R^b_n = c_b PiC(tau_n(H^b))` for every `n` and `b`. Source: DR s.3. Status: PROVED given
GI Lemma V and the Stage-1 definition `Pi_0 = PiC` through `Y^C` (MSG s.1). Depends on: GI Lemma V; D7.

**`prop:trm-dii` (prop) D(ii) (M3).** For all `b, k`: the `log p`-coefficient of `Phi_k^b` is `c_b Sum_i Pi_i(tau_k(H^b)) (-x)^i/i!`.
Hence for `s in S_w`: `Lambda(s)` is log-free iff `Pi_i(s) = 0` for all `i >= 0`, i.e. iff `s in K_0` (the `x^i` are independent);
and the log coefficient of the TR leg `-Lambda(s)` is `-Sum_i Pi_i(s)(-x)^i/i!` (SYN2 s.2 step 7.3).
* Source: statement DP s.1; proof of record DR s.3; SYN2 M3. Status: PROVED (DR s.3), given GI Lemma V and `Pi_0 = PiC` via `Y^C`;
  VERIFIED (critic Test A `k <= 100` over Q; DR D1-D3; `cr2_legdict` L3 `k <= 24`).
* Corrections: O6, O26 (the Stage-1 label "PROVED (s.6)" was a sketch: one sentence "LG dictionary"; proof of record is DR s.3);
  DR Correction 4(ii) (the closed form of `PiC` is an input).
* Depends on: `lem:trm-lemG`, `lem:trm-resid`, `Pi_i(tau_k H^b) = PiC(tau_(k-i) H^b)` (MSG s.1).

## s.2.2 The lattice `K_0 = L` (D(iii)); all from DP, refereed by DR ("no error, no gap")

**`lem:trm-count` (lem) Counting (DP Lemma 0).** For fixed `m` with `5m <= w` the weights `{1,2,3,4}` of `b = 0, 4, 3, 2` meet every class
mod 4 once; so `S_w` has `d(w)` triples with `b != 1` (the class-0 candidate `b = 2` fails exactly when `w = 5m`) and `e(w)` triples with
`b = 1`. Status: PROVED; VERIFIED (T1 `w <= 400`; DR R0 `w <= 4000`).

**`lem:trm-classL` (lem) The class L (DP Lemma 1).** (a) `Q[p,q] = Q[p,x] = Q[v][p] + Q[v][x]`, and `Q[p,x]` is free over `Q[v,x]` with
basis `1, p, p^2, p^3, p^4`. (b) `Q[v][p^(+-1)] cap Q(v)[x] = Q[v][x]`. (c) `L cap Q[v][p^(+-1)] = Q[p,q]`: a Laurent polynomial with
`Q[v]` coefficients is in `L` iff every term `v^m p^j` has `j >= -m`. (d) `L/Q(v)[x] = Q[p,q]/Q[v,x] = Q[v][p]/Q[v] = (+)_(i=1..4) Q[v,x] p^i`,
free over `Q[v]`, of dimension `d(w)` in weight `w`; byproduct `floor(w/4) + 1 = d(w) + e(w)`. Status: PROVED.
Correction: DR Correction 6 (the injectivity in (a) uses the algebraic independence of `p` and `x`: `v = px - p^5/4`; one line).

**`lem:trm-closure` (lem) Closure of M (DP Lemma 2).** (i) `I(p^a q^c) = -p^(a+4) q^c/(a-c+4) + p^a q^(c+1)/(a-c-1)`, a term with zero
denominator read as `-v^c log p`, resp. `+v^(c+1) log p`. (ii) For `g in Q[v,x]` with `x`-primitive `G`:
`I(log p . g(x)) = -log p . G(x) + Int G(x) dp/p`, where `Int p^a q^c dp/p = p^a q^c/(a-c)` (`a != c`), `= v^c log p` (`a = c`).
(iii) `I` is `Q[v]`-linear, `I(M) in M`, `I(Q[v,x]) in Q[v,x]`; so every `Phi_k^b in M`; moreover `M` is a `Q[v,x]`-module.
Status: PROVED. Correction: DR Correction 6 ("because `G(0) = 0`" unnecessary: the identity holds for every primitive; `g_0 != 0` occurs
when 5 divides the `x`-degree, e.g. `G = x^5/5`, `g_0 = v^4/4`).

**`lem:trm-modv` (lem) Reduction mod v (DP Lemma 3).** `Lambda (x) Q[v]/(v) : F/vF -> N/vN` is an isomorphism in every weight.
Status: PROVED (DR: bases, the relation `q = -p^4/4` mod `Q[v,x] + pq Q[p,q]`, closed images and triangularity checked).

**`thm:trm-lattice` (thm) D(iii) (M4).** For every `w >= 0`:
(a) `S_w` has `d(w)` triples with `b != 1` and `e(w)` triples with `b = 1`;
(b) for `s in K_0 cap S_w`, `Lambda(s) in Q[p,q] = Q[v][p] + Q[v][x]` (inside `L`), exactly, not only modulo functions of `x`;
(c) `Lambda : K_0 cap S_w -> (L/Q(v)[x])_w` is a bijection; both sides have dimension `d(w)`;
(d) `Lambda : S_w -> N_w` is a bijection; `N` is a free `Q[v]`-module with basis `{Lambda(tau_k(H^b))}`.
* Scope: coefficients in `Q[v]` are needed; over `Q(v)` the statement fails (`tau_0(H^2)/v` has leg `1/p`; DR K7: 42 of 222 log-free
  combinations outside `L` once `m = -1` is allowed).
* Source: DP s.1-s.5; DR s.0, s.2, s.6. Status: PROVED for all weights (no TR, no quoted theorem); (a), (d) need nothing; (b), (c) need
  D(ii). VERIFIED exact `w <= 400`, two primes `w <= 1000`.
* Corrections: O1/DR Corr. 1 (the `c_b` are not cosmetic for `K_0` as a subspace: with raw legs the log-free kernel differs from `K_0`
  at 390 of 401 weights; the "16200 `K_0` vectors" are raw-leg kernel vectors); O2/DR Corr. 2 (`dim K_0 cap S_w = d(w)` comes from (d)
  or Cor. 4 together with D(ii), not from (a)); K11 (D(iii) is PROVED; LL's "VERIFIED to weight 240" is outdated).
* Depends on: `lem:trm-count`, `lem:trm-classL`, `lem:trm-closure`, `lem:trm-modv`, `prop:trm-dii` (for (b), (c)); graded Nakayama.

**`cor:trm-split` (cor) DP Cor. 4.** `S_w = (K_0 cap S_w) (+) span{v^m tau_k(H) : 5m + 4k = w}`. Status: PROVED; uses D(ii) as stated
with `K_0` (K3, O3, DR Corr. 3; D(ii)-free only with "log-free combinations" in place of `K_0`). Every `s in S_w` is
`s_0 + Sum a_(m,k) v^m tau_k(H)` with `s_0 in K_0`, the `a_(m,k)` read off from the log component.

**`cor:trm-Ksplit` (cor) DP Cor. 5.** `K cap S_w = (K_0 cap S_w) (+) Q v^(w/5) tau_0(H)` if `5 | w`, and `K cap S_w = K_0 cap S_w`
otherwise. So on `K` the only legs beyond `Q[p,q]` are multiples of `log p`. Status: PROVED given D(ii) (DP s.9; DR s.2 checked:
e.g. `K cap S_5 = S_5` contains `tau_1(1)`; `K_0 cap S_9 = <v tau_0(H^2), tau_2(1) + v tau_1(H)>`).

**`prop:trm-coeff` (prop) Coefficient ring (M5).** Over `Q[[u]]`, `K_0` is spanned by the frames `frame(k, a)` (`a != 1`) with
`f_j = const . v^((4(k-j)+c_a+1)/5)`, exponent `>= 1`; so every insertion of quintic Theorem R is a `Q[[u]]`-combination of elements of
`K_0 cap S_w` and `tau_0(H)`. With DR Correction 5(iii): the generators `kappa(X, k)` have `Q[v]` coefficients, and the `v^m kappa(X,k)`
of weight `w` form a basis of `K_0 cap S_w` (independent, in `K_0 cap S_w`, `d(w)` of them).
* Source: SYN2 M5; DR s.2 Correction 5. Status: PROVED (theory Lemma 0 plus the exponent count); VERIFIED (`cr2_legdict` L4, `k <= 12`).
* Depends on: RD Lemma 0; `thm:trm-lattice`.

**`rmk:trm-h2leg` (rmk) Sensitivity to the primary legs.** Replacing `Phi_0^2 = q` by `x` kills injectivity at `w = 4`; for
`Phi_0^2 = a p^4 + b q`, Lemma 3 predicts (d) iff `a != b/4` (i.e. iff `Phi_0^2` is not a multiple of `x`), VERIFIED for 7 choices,
`w <= 40`; by weight, `Phi_0^2` is the only primary leg not fixed up to scale (`M_4 = <p^4, q> (+) log p . <x>`). Source: DR s.2
Correction 4(i); SYN2b O4. Status: remark (prediction from Lemma 3 + VERIFIED).

## s.2.3 L-regularity, the Gap Theorem, log legs

**`thm:mlk-lreg` (thm; quoted) L-regularity (E9, M7).** On the family `x = z^4 + c3 z^3 + c2 z^2 + c1 z + V/z`, `y = z`, for `2g - 2 + n > 0`
every correlator whose `n >= 1` legs all lie in `Q[c,V][z] + Q[c,V][x]` is a polynomial in the parameters (so regular at `V = 0`).
On the line (form D1): a pairing of `omega_{h,n}` with all legs in `Q[v][p]` lies in `Q[v]`.
* Source: XY s.0 item 4, s.3 (ML1); XYR s.3; S01 E9; SYN2 M7. Status: PROVED mod arXiv:2212.00320 Thm 7.1 (statement and hypotheses
  checked by XYR s.2, plus an off-line test); VERIFIED in every falsifier-(b) pairing, `h <= 6`.
* Depends on: ABDKS Thm 7.1 (7.9); the kernel sign `omega^ABDKS_{g,n} = (-1)^n omega^engine_{g,n}` (XYR s.2); `lem:trm-xfun`; the scaling
  between the forms of D1 (s.7 item W5).

**`thm:mlk-gap` (thm; quoted) Gap Theorem (M8).** For `h >= 2`, in the engine convention, `F_h^TR = B_2h/(2h(2h-2)) v^(2-2h)` on `x = W`
(`F_2^TR = -1/240`); hence frozen `F_h = (-5)^(1-h) F_h^TR = c_h v^(2-2h)`. General form (TRM s.4.1, "Theorem G" there): for
`x = P(z) + v/(z - a)`, `P` of degree `r >= 1`, `v != 0`, simple ramification, `y = z`: `F_g = B_2g/(2g(2g-2)) v^(2-2g)`, `g >= 2`.
* Source: SYN2 M8; TRM s.0 item 3, s.4.1; XY s.2.3; XYR s.2. Status: PROVED by two routes: (1) arXiv:2502.18115 Cor. 1.2, formula (1.3)
  (TRM s.4.1 proof; its printed example s.3.2 has a sign slip); (2) arXiv:2212.00320 Thm 7.1 + a Borel sum (XY s.2.3; XYR s.2 short
  proof). VERIFIED `h <= 6`.
* Corrections: LR F5, SYN2b O13 (state the convention of the quoted corollary: `F_h` here is `-F_EO^print`; TRM s.1 says Hock's kernel
  and `F_g` (1.1) are "identical to ours"); TRM's "Theorem G" renamed "Gap Theorem" (S01 C9).

**`lem:mlk-vanper` (lem; quoted) Vanishing period.** For `g >= 1`: `Int_oo^0 omega^engine_{g,1} = (B_2g/2g) V^(1-2g)` on the family of
`thm:mlk-lreg`; hence `<-log z>_g = -Int_oo^0 omega_{g,1} = -(B_2g/2g) V^(1-2g)` (orientation: XYR s.2 "F_1" item). Source: XY s.2.3;
XYR s.2. Status: PROVED mod arXiv:2212.00320 Thm 7.1; VERIFIED `g = 1, 2, 3` at generic points. Not used by the sources' chain; it is a
written alternative for the base `h = 1` of Theorem V(ii) (s.7 item G1).

**`lem:mlk-fixedx` (lem) Fixed p versus fixed x (LL s.2.1).** At fixed `x`, `d_v p|_x = -p/(p^5 - v)`; for any family of forms
`d_v|_x = d_v|_p - Lie_xi` in each variable; `d_v(y dx)|_x = -dp/p = Int_{q in gamma} B(., q)`: a third-kind variation (residue `-v` of
`y dx` at `p = 0`, `+v` at `oo`), "`Lambda = 1` on `gamma`" in the engine convention. Status: PROVED (elementary); LR re-derived.

**`lem:mlk-ibp` (lem) Integration by parts (LL s.2.2).** (a) L-slot: `oint Phi Lie_xi w = -oint (xi . dPhi) w = -<(xi dPhi)_+, w>`.
(b) Log slot: `Int_gamma Lie_xi w = [i_xi w]_0^oo = 0` for `w` holomorphic at 0 and `O(dp/p^2)` at oo. Status: PROVED; LR checked.

**`thm:mlk-var` (thm; quoted) The variational formula (Q1), in the form used.** For the family D1 (fixed `B`, variation at fixed `x`
with `d_v(y dx)|_x = Int_gamma B(., q)`), for `2h - 2 + n >= 0`, `(h,n) != (0,0)`:
`d_v|_x omega_{h,n}(J) = Int_gamma omega_{h,n+1}(J, q)` and, for `h >= 2`, `d_v F_h = Int_gamma omega_{h,1}`.
* Source: LL s.2.3; LR s.2 (dictionary). Status: quoted: arXiv:math-ph/0702045 Thm 5.1, (5-14) for `g + k > 1` and (5-15) for `g >= 1`,
  with the temperature variation of its s.5.3, (5-27)-(5-30), read through the dictionary `omega_{h,n} = (-1)^n W^print_{h,n}`,
  `F_h = -F^print_h`, `Lambda^print = -1` (LR s.2). VERIFIED as a tensor identity (14 cases worker, 7 referee). The case `(0,1)` is the
  hypothesis; the case `(0,2)` is re-proved directly (`lem:mlk-genus0`).
* Corrections: F1/O7 (LL's "our kernel is EO's" and "s.5.1" are wrong: the kernel is opposite to the printed (4-3); the temperature
  variation is s.5.3; the "form used" is correct only through the dictionary; the two sign discrepancies cancel).
* Hypotheses as checked (LL s.2.3; LR s.2): simple branch points with `dy != 0` for all `v != 0`; variation at fixed `x`, `B` the
  Bergman kernel of the fixed P^1; `Omega/dy` regular at the branch points (EO (5-4)); the path `[0, oo]` joins the two poles of `y dx`
  and avoids the branch points (EO (5-13), footnote 8: no hard edges); the integrals converge (`omega_{h,n+1}` holomorphic at 0 and oo);
  residues at the moving branch points commute with `d_v|_x`.

**`lem:mlk-genus0` (lem) Genus zero (LL s.2.4).** `d_v|_x B = Sum_b r_b e_{b,1}(p1) e_{b,1}(p2)`, `r_b = 1/(5 c_b^3)`, `e_{b,k} = dp/(p - c_b)^(k+1)`;
`omega_{0,3} = -Sum_b e_{b,1} e_{b,1} e_{b,1}/x''(c_b)`; `Int_gamma e_{b,1} = -1/c_b`; hence `Int_gamma omega_{0,3}(p1, p2, .) = d_v|_x B`,
the case `(0,2)` of (Q1) with the engine sign. Explicit (0,3) pairings (`Sum_b c_b^N = 5 v^(N/5)` if `5 | N`, else 0):
`<Phi1, Phi2, -log p>_0 = (1/5) Sum_b c_b^-3 Phi1'(c_b) Phi2'(c_b)`; `<Phi1, (-log p)^2>_0 = -(1/5) Sum_b c_b^-4 Phi1'(c_b)`;
`<(-log p)^3>_0 = 1/v`; `<Phi1, Phi2, Phi3>_0 = -(1/5) Sum_b c_b^-2 prod Phi_i'(c_b)`. Every (0,3) pairing with an L-leg lies in `Q[v]`.
Status: PROVED; LR re-derived. Correction: F1/O13 (this sign of `omega_{0,3}` is opposite to EO Thm 4.1 (4-20); say so).
