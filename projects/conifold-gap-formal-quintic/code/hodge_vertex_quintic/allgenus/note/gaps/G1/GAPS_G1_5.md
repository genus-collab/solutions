# Gap cluster G1, part 5 of 5: item 5 (the model), consequences for the design, reading record

Continuation of `GAPS_G1_4.md`. Model conventions of OUTLINE_3 s02 items 10 (`def:set-model`, `prop:set-modrows`, `def:set-M0`,
`def:set-modcohft`): `W_v(p) = p^4/4 + v/p`; rows `Y_b = eps_b I_(nu_b)`, `nu = (0,-1,3,2,1)`, `eps = (1,-1,-1,1,1)`; model system
`z dY_b/dv = bm_b Y_(b+1)`, `bm = (1, -1/v, 1, -1, -1)`; `varsigma = v^(1/5)`, `t = z/varsigma^4`; colour-0 solution `Y_b = e^(U/z) varsigma^(nu_b - 1)
Fs_b(t)`, `U = -(5/4) varsigma^4`, `Fs_0(0) = 1`, `Fs_b(0) = eps_b` (`prop:set-modrows` (b), UR s.2 Prop. S); `Rhat_b^(i)(t) := zeta^(i(b+1))
Fs_b(zeta^(-i) t)`; `Psi_0 := Rhat(0)`; `M_0(t) := Psi_0^(-1) Rhat(t)`; `eta^(ab) = (1/5) delta_(a+b = 3 mod 5)`.

## 5. Item 5

**Search result (FOUND only as statements).** Symplecticity of `M_0`: LG_MODEL s.3 "Its output is symplectic through `zhat^40`" (VERIFIED);
UR s.7 table "symplectic, Z/5-covariant" (asserted); MSG2 "Symplecticity is automatic for the actual `R' = R M^-1` (both factors are
symplectic R-matrices)" (asserted). Homogeneity of the model CohFT: MSG1 s.1 "Everything is homogeneous ... a genus-h correlator is
`const * v^e`" (stated); RD1 s.1.2 "(PROVED, LG_MODEL s.3, UNIFORM_REDUCTION s.1)"; but UR s.1 and LG_MODEL s.3 prove only the homogeneity
of the ROWS (`x -> kx, v -> k^5 v, z -> k^4 z`), not of the CohFT. Both are therefore PROVED here.

### 5.1 Lemma (symplecticity of the model R-matrix): PROVED

*Statement.* `Sum_(a,b) Rhat_a^(i)(-t) eta^(ab) Rhat_b^(j)(t) = delta_ij`; equivalently `M_0(-t)^T M_0(t) = 1`, i.e. `M_0(t)^(-1) = M_0(-t)^T`.

*Proof.* (1) *Colour-`i` columns.* In the variable `varsigma` the system reads `(z/(5 varsigma^4)) dY_b/d varsigma = bm_b(varsigma^5) Y_(b+1)`,
which is invariant under `varsigma -> zeta^(-i) varsigma`. Applying this to the colour-0 solution, and using `nu_b = -b mod 5` (so
`zeta^(-i(nu_b - 1)) = zeta^(i(b+1))`), `Y^(i)_b := e^(U_i/z) varsigma^(nu_b - 1) Rhat_b^(i)(t)`, `U_i = -(5/4) zeta^i varsigma^4`, is a formal solution
for each `i`. (2) *Self-adjointness.* With `B := Sum_b bm_b E_(b,b+1)`, `(eta^(-1) B)_(ac) = (1/5) bm_(3-a) delta_(c, 4-a)` is symmetric iff
`bm_(a-1) = bm_(3-a)` for all `a` in Z/5, which holds for `bm = (1, -1/v, 1, -1, -1)` (`bm_4 = bm_3`, `bm_0 = bm_2`). (3) *Constancy.* For
`P_ij := Sum_(a,b) Y^(i)_a(v, -z) eta^(ab) Y^(j)_b(v, z)` (the first factor solves the system with `-z`), `d P_ij/dv = z^(-1) Y^(i)(-z)^T
(eta^(-1) B - B^T eta^(-1)) Y^(j)(z) = 0`. (4) *Reduction.* On `eta`-pairs `nu_a + nu_b = 2`, so the powers of `varsigma` cancel and `P_ij =
e^((U_j - U_i)/z) S_ij(t)`, `S_ij(t) := Sum_(a,b) Rhat_a^(i)(-t) eta^(ab) Rhat_b^(j)(t)` in `Q(zeta)[[t]]`. Using `dU_k/dv = -zeta^k varsigma^(-1)`
and `dt/dv = -(4/5) t/v`, `dP_ij/dv = 0` is equivalent to `(zeta^j - zeta^i) S_ij(t) + (4/5) t^2 S_ij'(t) = 0`. For `i != j` the lowest
coefficient of `S_ij` would have to vanish, so `S_ij = 0`; for `i = j`, `S_ii' = 0`, so `S_ii = S_ii(0) = (Psi_0^T eta^(-1) Psi_0)_ii = 1`
(`def:set-M0`: `(1/5) Sum_a eps_a eps_(3-a) zeta^(i(a+1) + j(4-a)) = delta_ij`). (5) *Teleman form.* `Psi_0^T eta^(-1) Psi_0 = 1` gives `Psi_0
Psi_0^T = eta`, so `M_0(-t)^T M_0(t) = Rhat(-t)^T (Psi_0 Psi_0^T)^(-1) Rhat(t) = Rhat(-t)^T eta^(-1) Rhat(t) = 1`. QED

Inputs: only `prop:set-modrows` (a)-(b) (the model system and its WKB form). Formal throughout (no integrals): `Y^(i)` are formal WKB
columns, `d/dv` acts by `d e^(U/z)/dv = (U'/z) e^(U/z)` and the chain rule in `varsigma`. Corollary: the outer `R' = R M^(-1)` is symplectic
(product of symplectic series), the assertion of MSG2.

### 5.2 Lemma (homogeneity of the model CohFT): PROVED

*Definitions used* (`def:set-modcohft`). Colour space `C^5` with orthonormal basis `e~_i`; TFT `omega^m_v(e~_(i_1), .., e~_(i_n))_h =
delta_(i_1 = .. = i_n = i) (Deltahat_i^(1/2))^(2h-2+n)`, `Deltahat_i^(1/2) = zeta^(4i) varsigma`; `Omega' := M_0(z/varsigma^4).omega^m_v` (PPZ15 Def. 2.13);
`Psi^m_v : (C^5, delta) -> (H, eta)` with lowered entries `eta(phi_b, Psi^m_v e~_i) = varsigma^(nu_b - 1) eps_b zeta^(i(b+1))` (`= Q^m(v; z = 0)`),
an isometry (the computation of s.5.1 (4) with `nu_a + nu_b = 2`); `Omega^mod_v := (Psi^m_v)_* Omega'`, `(A_* Omega)_(g,n) = Omega_(g,n) o (A^(-1))^((x)n)`.
(U) *Unit:* `(Psi^m_v)^(-1) phi_0 = varsigma^(-1) Sum_i zeta^i e~_i = Sum_i Deltahat_i^(-1/2) e~_i`, the unit of `omega^m_v`; so `Omega^mod_v` has unit
`phi_0` (PPZ15 Prop. 2.12; used in GAPS_G1_4 s.4 (ii)). *Scope:* the value `Deltahat_i^(1/2) = zeta^(4i) varsigma` is taken from
`def:set-modcohft` as written in OUTLINE_3 (S15 (I4)); if the note instead defines `omega^m_v` through the Hurwitz Frobenius manifold of
DNOPS (8-1), this value is the consistency statement S15 (I4), which the lemma below uses and does not re-derive.

*Statement.* (a) For every cohomological degree `D`, `Omega^mod_(v; h,n)(phi_(b_1), .., phi_(b_n))_[D] = varsigma^(2h-2+Sum_j nu_(b_j) - 4D)
Omega^mod_(1; h,n)(phi_(b_1), .., phi_(b_n))_[D]` (`varsigma = 1` on the right). (b) Hence `<prod_j psi^(k_j) phi_(b_j)>_h(v) = C v^e`, `C =
<..>_h(1)`, `e = 2 - 2h + Sum_j wt(k_j, b_j)`, `wt(k, b) = (4k + nu_b - 4)/5`, with `v^e := varsigma^(5e)`. (c) `C = 0` unless `e` is an integer.

*Proof.* (a) First for `Omega'` in the basis `e~_i`, graph by graph in PPZ15 Def. 2.2 with Def. 2.13 (`R.omega = R T omega`, `T(z) = z 1 -
z R^(-1)(z) 1`, `R(z) = M_0(z/varsigma^4)`). Every ingredient is `varsigma`-free except: a leg monomial `psi^k` of `R^(-1)(psi)` carries `varsigma^(-4k)`;
an edge `(1 - R^(-1)(psi') R^(-1)(psi'')^T)/(psi' + psi'')` is `varsigma^(-4)` times a series in `psi'/varsigma^4, psi''/varsigma^4`, so its monomial
`psi'^a psi''^b` (degree `a + b + 1` after the boundary push-forward) carries `varsigma^(-4(a+b+1))`; a kappa-leg `T(psi)` has its `psi^m`
coefficient proportional to `varsigma^(-1) varsigma^(-4(m-1))` (`1 = varsigma^(-1) Sum_i zeta^(-4i) e~_i`) and contributes degree `m - 1` after push-
forward, but also raises the valence of its vertex by one; a vertex of genus `h_v` and total valence `n_v + m_v` (`m_v` kappa-legs) carries
`varsigma^(2h_v - 2 + n_v + m_v)`. Since `Sum_v (2h_v - 2 + n_v) = 2h - 2 + n` and `1 + (-1 - 4(m - 1)) = -4(m - 1)` for each kappa-leg, every term
of degree `D` carries `varsigma^(2h - 2 + n - 4D)`. Then `(Psi^m_v)^(-1) phi_b = varsigma^(nu_b - 1) Sum_i eps_b zeta^(i(b+1)) e~_i` contributes
`varsigma^(Sum_j (nu_(b_j) - 1))`. (b) Integrating against `prod psi_j^(k_j)` selects `D = 3h - 3 + n - Sum k_j`; the exponent becomes
`2h - 2 + Sum nu_(b_j) - 4(3h - 3 + n - Sum k_j) = 5(2 - 2h + Sum_j (4k_j + nu_(b_j) - 4)/5) = 5e`. (c) The flat R-matrix is `Psi^m_v M_0(z/varsigma^4)
(Psi^m_v)^(-1) = Q^m(v; z) Q^m(v; 0)^(-1)` (S15 Lemma 1: frame changes conjugate), where `Q^m(v; z)` is the matrix of the columns `e^(-U_i/z)
Y^(i)`. Replacing `varsigma` by `zeta^(-1) varsigma` permutes these columns (`Y^(i) -> Y^(i+1)`, s.5.1 (1)) and the TFT data (`Deltahat_i^(1/2) ->
Deltahat_(i+1)^(1/2)`), so `Omega^mod_v` in the flat basis does not depend on the choice of the fifth root; by (b), `(zeta^(-5e) - 1) C = 0`. QED

## 6. Consequences for the design (for the editor; nothing here edits a design file)

1. *s02 item 2.* Replace (S1) by (M) (s.1.1 wording; pin s.1.2-s.1.4). Add `prop:set-cyclicLR` (register in LABELS) with the proof of
   s.2.8 and Lemmas A-F as local lemmas (prefix `set`; about 1.5 pp.), and a corollary (proposed `cor:set-cyclicPF5`) = Corollary 2.10
   (a)-(b): `Delta_i`, `D u_i`, and "`C_0 v_0` solves PF5 for every solution of the QDE". The `\draftnote` of R1 fix item 1 can go.
2. *`def:set-rows`.* Define `Phi_i := C_0 Q^(i)_0/(zeta^i L)` and `A_n := [z^n](C_0 Q^(0)_0)`; the WKB/PF5 property becomes Cor. 2.10 (b).
3. *`thm:set-qrr`.* Delete `\draftnote` (a) (s.3.1). Replace (b) by Steps 1-5 of s.3.2 with Teleman v3 Theorem (2) in (1.15), 8.7
   Proposition (i)-(ii), 8.10 Proposition, PPZ15 Prop. 2.14, LP18 (26), GJR18 Thm 4.1; keep ONE `\draftnote`: "injectivity of the
   classification read from Teleman 8.7 (ii); the rank-one elementary route uses Madsen-Weiss, not re-read".
4. *`lem:set-dilaton`.* Proof of s.4 with PPZ15 Def. 0.3 (iii), Prop. 2.12, Def. 2.13 and Arbarello-Cornalba (1.7); state the hypothesis
   `kappa` in `1 + m`; ask s04/s07 to state `d(0) = 1` where it is applied. The `\draftnote` "no pinned citation" can go.
5. *`def:set-M0`, `def:set-modcohft`.* Add two lemmas (proposed `lem:set-M0symp`, `lem:set-modhom`) with the proofs of s.5.1, s.5.2.
6. *R7 table.* Remove `thm:set-qrr` (a)-(b), `prop:set-cyclicLR`, `lem:set-dilaton`, `def:set-M0` symplecticity, `def:set-modcohft`
   homogeneity from the list of open steps; the conditionality of the main theorem then names (M) with its pin; (Q) and (QRR) unchanged.
7. *BIBPLAN.* MR1408320 (Givental 1996) is now used (Cor. 11.7, Thm 11.6). New keys: Arbarello-Cornalba alg-geom/9406008v1 (journal
   reference to be confirmed: the arXiv page lists none); optionally Guo-Ross 1703.06955v1 (Compositio Math. 155 (2019) 995-1024, from the
   arXiv page) for the literature remark of s.2.11. Do not add Zagier-Zinger as a proof citation (not read; not needed).
8. *Notation (R8).* `d_T` (derivative in the mirror coordinate `T`) and the formal `log q` (here `ell`) need local letters: `ell` and `d_t`
   are taken in s03 (NOTATION table 1).

## 7. Reading record (what was read in full text, as images of the arXiv PDFs)

* Givental alg-geom/9603021v2: p.1-3, p.22-29 (s.8 end, s.9 set-up, Thm 9.1, Cor. 9.2-9.4, Thm 9.5, Prop. 9.6-9.8), p.35-42 (s.11: Main
  Lemma, Cor. 11.1, Props 11.2-11.5, Thm 11.6, Cors 11.7-11.8 and remarks). Not read: s.1-s.8 except as named, p.30-34, p.43-48.
* Coates-Givental math/0110142v2: p.1-9 (intro, s.1-s.7 up to Thm 2). Not read: proofs of Thm 1-2, s.8, appendices.
* LP18 (1803.01409v2): p.1-8, p.10-19, p.21-23, p.25-26, p.38-39 (references). GJR18 (1812.11908v1): p.1-2, p.13-14, p.18-19, p.26-28, p.45-46.
* Guo-Ross 1703.06955v1: p.1-5, p.15-21. Teleman 0712.0160v3: p.1-2, p.4-11, p.39-41. PPZ15 1303.1043v3: p.6-10, p.19-20, p.24-30.
  Arbarello-Cornalba alg-geom/9406008v1: p.4-6.
* NOT read: Lian-Liu-Yau; Zagier-Zinger 0710.0889; Kim-Lho; Ciocan-Fontanine-Kim; Faber-Pandharipande (Mumford's relation used as
  recorded in LIT s.1.4); Madsen-Weiss; Witten (for `Int_(Mbar_(1,1)) psi = 1/24`); Guo-Ross 1611.08876; Pandharipande 1712.02528.
* Project sources searched: FQ2 s.1, s.5; fq_rows.py, fq_genus0.py, genus0_check.py docstrings; PLAN_2026-11-03 "Key facts"; LIT s.0-s.1;
  UR s.1-s.2 and s.7 table; LG_MODEL s.0-s.4 headers and statements; MSG1 s.0-s.1; MSG2 (grep); RD1 s.1.2; S15 s.1-s.2; PS s.4.
  Closest project text for item 2: the docstring of `fq_genus0.py` ("reading off the z^0 term of z d_T S(H^a) = S(H*H^a) gives the
  structure constants, and at the zeta point Z/5-covariance forces H*H^a = kappa_a H^(a+1 mod 5) (checked)"), i.e. the Birkhoff
  recursion and the grading as a computational recipe, with the values checked to `q^12` (genus0_check.py; FQ2 s.1, s.5). No proof.
