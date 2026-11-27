# GAPS_G2 part 6: evidence (s.10), corrections and notes for the editor (s.11), dependencies (s.12), files (s.13)

Continuation of `GAPS_G2_5.md`.

## 10. Evidence (VERIFIED; exact; none of it is used by a proof)

All scripts are stdlib Python 3.14, run as `py -I -B <script> <args>` from this folder (each adds its own folder to `sys.path`;
they import only `g2_core.py` and each other, no project module). Logs print no paths. Every log ends with `RESULT: ALL CHECKS PASSED`.

| script / log (args) | what is tested | range, field | failing controls (each fails as it must) |
|---|---|---|---|
| `g2_linalg.py` / `g2_linalg.log` (6) | s.3 (a)-(e) and the lifting formula of s.4.1 at random points `b~` in `(Q^x)^5`: `det(x - A~) = x^5 - Bprod`; `rank ad A~ = 20`; commutant = `span(1, .., A~^4)`, `A~^5 = Bprod 1`; Gram determinant `5^5 Bprod^4`; `X - pr(X)` in `im ad A~`; the `calY`-form; disjoint supports of the `calY^j`; `[A~, calZ] = X - pr X` | 6 points, exact `Q` (54 checks) | factor `1/4` in place of `1/5`: `X - pr X` not in the image (rank 21); `A~^2` not in the image |
| `g2_lemma2.py` / `g2_lemma2.log` (14) | s.4.2 (a): `pr((calY^j)') = -(j/5)(Bprod'/Bprod) calY^j`, `pr(dlog Theta_U calY^j) = (tr/5) calY^j`, `j = 0..4`, with `pr` computed from the trace form by matrix products and traces (not by the closed formula); residue `j/5`; pole bound of s.3 (f) on a holomorphic test matrix | quintic and pulled-back model, frame (0,0), mod `u^13`, exact `Q` | `j/4` in place of `j/5`; a non-diagonal matrix in place of `dlog Theta_U` |
| `g2_gauge.py` / `g2_gauge.log` (15 32) | the construction of s.5.3 in graded form (`lem:gau-lift` with `rho_1 = 0`; `lem:gau-scal` with `gamma = 1` at `k = 0`, `0` at `k = 5q > 0`); the log obstruction matrix of PF4 recomputed (`gamma = 7/10`); (EP), (N0); gauge equation in every entry and known coefficient; holomorphy; `P_0` diagonal units; `P_11(0,z) = 1`; G8 (`p_(0,a)(0) = 1`, `prod p_(0,a) = 1` to `u^30`); `p_(0,0) = tau^(-1/4)/varpi_0` | frames `(0,0)`, `(1/3,-2/7)`; `z^0..z^15`; `u`-precision 31 down to 19; exact `Q` | `2 P_can`: product 32; deformed model `b~_0 -> 2 b~_0`, `b~_1 -> b~_1/2`: product `1/16` (as s.7 predicts) |
| `g2_G.py` / `g2_G.log` (15 32 10) | (G) for the constructed `P_can` in every slot `z^-13..z^3`, mod `u^10` (85 slots per frame); `[z^-1] Y^C_Q = -t_c e_0`; `P_can c(z)` gives `c(z) Y^C_Q`; the canonical values of the sources: `pi(z) = (127/3150) z^3 + (1416472237/725359556250000) z^8` (frame (0,0)), `-(239/450) z^3 - (20846188763/725359556250000) z^8` (frame `(1/3,-2/7)`), `p_(1,b)(0) = (8/45, 0, 8/45, -1/15, -1/15)` resp. `(-7/45, 0, -7/45, -1/15, -1/15)` (GI s.5; VT `v2_run.log`) | both frames, exact `Q` | `P_can (1 + z^5/7)`: `[z^5 u^0]` slot 1 is `-1/7`, so `P Y^C != Y^C_Q` (the uniqueness of G9) |
| `g2_unique.py` / `g2_unique.log` (16 24 2) | the uniqueness clause INCLUDING non-graded gauges, by a route independent of the construction: all equations of `GE(A, A~)` that involve only the unknowns `[u^m](P_k)_(ab)`, `k <= 16`, `m <= 24`, split by the Z/5-degree; the projection of the solution space onto the triangle `k <= 15`, `m <= 22 - k` has dimension `#{k <= 15 : k = delta mod 5}` = `4, 3, 3, 3, 3` for `delta = 0..4`, and `z^delta P_can` solves every equation (so the projection is spanned by `P_can z^r Q[[z^5]]`) | both frames, mod `2^61 - 1` | unequal products (`b~_3 -> (1+u) b~_3`): dimension 0 on the box `k <= 5`, `m <= 10` (equal products on the same box: `2, 1, 1, 1, 1`); residue equations dropped: 29 and 28 (not 4, 3); windows touching the truncation (margins 0, 1): 6 and 4 (not 3), i.e. the method does see spurious solutions |
| `g2_colours.py` / `g2_colours.log` (10) | `lem:gau-colmod`: the T-system for colours 0..4, `F^(i)_b(t) = zeta^(ib) Fs_b(zeta^-i t)`, pivot `4n zeta^(4i)`, `Fs_0 = 1 - (3/20)t + (9/800)t^2 + (11571/80000)t^3 - (2415609/6400000)t^4 + ..`, `Fs_b(0) = (1,-1,-1,1,1)`; `w_i = zeta^i w_0` (`lem:set-hrot`); gradedness of the stored `P_can` entries (structural) | to `t^10`, exact `Q(zeta)` | normalised rows `zeta^(i(b+1)) Fs_b(zeta^-i t)` are not the T-system solution; `{zeta^i - zeta^j} U {-5}` not covariant; `P_can (1 + z)` not covariant (structural) |

Remarks on the evidence. (1) With the choices of `g2_gauge.py` (`rho_1 = 0` in the lift, `gamma = 0` at `k = 5q > 0`) the
construction gives `P_11(0, z) = 1` exactly, i.e. it produces `P_can` directly under (N0): `(calZ_(5q))_(11) = rho_1 = 0` and
`(P_0)_(11)(0) = d^U_1(0) = 1`. (2) The uniqueness test needs a window that stays two `u`-orders per `z`-order away from the truncation;
for unequal products the truncation artefacts travel about five `u`-orders per `z`-order (the inverse of `X -> A X - X A~` then has
poles in `u`), hence the separate box for that control. (3) Not computed here: the actual WKB columns of colours `i != 0` and the
series `D_i`; the covariance (C2) of the actual columns is VERIFIED in the sources only (FQ2 s.1; CS s.1).

## 11. Corrections and notes for the editor (author-facing; for `DRAFT_NOTES.md` / `SOURCE_ERRATA.md` as appropriate)

1. **`thm:gau-U`.** State 5.1 (a)-(d) (and (e) as a remark under `prop:gau-equalprod`). The `\draftnote`s G1, G2, G3, G4, G5, G7 of
   OUTLINE_4 can be replaced by the proofs of s.3-s.6; the draftnote "the construction is written over `C[[u]]` ... over `Q` is not
   discussed" by 5.5. Add `lem:gau-lead` before the proof of the uniqueness clause and `lem:gau-degree` before `prop:gau-graded`.
2. **`cor:gau-colours`.** The three-line computation (9.2) needs (C1), (C2), (C3). (C2) rests on `lem:gau-colrows`, whose input (i)
   (cyclic shape of the actual large-radius rows) is part of the open item R1 (`prop:set-cyclicLR`, cluster G1) and whose input (iii)
   contains (M) (cluster G1, D2) and the colour-`i` form of (QRR) (`thm:set-qrr` draftnote (a)). Recommended placement:
   `lem:gau-colmod`, `lem:gau-colrows` in `sec:gau-R0` after `lem:gau-wkb`; `def:set-rows` then says "covariance: `lem:gau-colrows`"
   instead of citing it as a fact. Until R1 is closed, `cor:gau-colours` carries one `\draftnote` naming `prop:set-cyclicLR`.
3. **`thm:gau-R0` (iv), DB D14, PS s.4.1, DD1 D1.8.** "`D_i` in `Q[[z]]`" holds for `i = 0`; for `i != 0`, `D_i(z) = D_0(zeta^-i z)` lies
   in `Q(zeta)[[z]]`, and the comparison ring of colour `i` is `Q(zeta)[L, 1/L][[z]]` inside `Q(zeta)((s))[[zhat]]`. Harmless after
   `D_0 = 1`, but the statement must say it.
4. **Integration constants of the `u_i`.** Not needed: `lem:gau-colrows` and 9.2 work with amplitudes, which see only `D u_i =
   zeta^i L`; the exponent remark of 9.2 uses only `u_i - u_i^c = zeta^i (u_0 - u_0^c)` (`lem:set-canon`). (A first draft of this note
   asked for a normalisation `u_i = zeta^i u_0`; that request is withdrawn.) What must be kept is PS-C8: `e^(u_i^c/z)` is not absorbed
   into `D`.
5. **`lem:gau-diag`.** "`= 1` for `P_can`" holds under (N0), i.e. in the programme family (`t_c'(0) = varpi_0(0) = 1`, `varpi_1 = u +
   O(u^2)`, `v = u + O(u^2)`); state the general value `p_(0,1)(0)^5 prod_a rho^U_a` (s.7) or restrict the statement explicitly. The
   second proof of s.7 (from `lem:gau-lead`) is shorter than GI s.5 steps 1-3 and gives `p_(0,0) = (det Theta_U)^(-1/5)` directly.
6. **`prop:gau-symp` step 3.** Cite 5.1 (b) for the pair `(A~, A~)` in the form "every solution of `z T' = A~ T - T A~` with integral
   coefficients is `t(z) 1`" (s.5.4, consequence 2); invertibility of `T_0` is then not needed.
7. **`def:gau-onepole` / `thm:gau-U`.** "`P_0` diagonal with unit entries" holds for every regular gauge (`P_0 = c(0) frak d Theta_U`),
   not only for the constructed one; "regular gauge" can be defined by `P_0 != 0`.
8. **`prop:gau-graded` step 1.** Use `lem:gau-degree` (no `zeta`, no averaging, works over `Q`); the `S`-conjugation and the average
   `(1/5) Sum_i Rot_i P = P^[0]` are equivalent formulations (6.4 (i)); keep at most one sentence on them.
9. **`rmk:gau-Uremarks`.** (a) now has its proof (5.6). (b) The spectral-curve reading stays a remark; its first half (`F[calY]` is a
   field, totally ramified) is s.3 (g) (Eisenstein), the rest is not used. "Regular semisimple" (UR s.3 (ii)) is true (s.3 (g)) but not
   needed: say so, or drop the words.
10. **New labels** (register in `LABELS.md`): `lem:gau-lead` (leading term of a solution of the gauge equation), `lem:gau-degree`
    (Z/5-degree of the gauge equation), `lem:gau-colmod` (colour covariance of the model columns), `lem:gau-colrows` (colour covariance
    of the actual columns). **New local symbols** (top of `s03`, local to `sec:gau-U` unless stated): `calB_a^(m)` (partial products;
    `\mathcal B^{(m)}_a`, extending `Bprod`), `calY := A~^-1` (`\mathcal Y`), `dlog Theta_U := Theta_U^-1 Theta_U'`, `frak d`
    (`\mathfrak d`), `pr`, `lambda_m` (local to the proof of `lem:gau-comm`), `kappa^U_j` (or another letter: `kappa_a` is taken),
    `P^\natural` (the constructed gauge), the Z/5-degree `P^[delta]`; for `sec:gau-R0`: `Rot_i`, `T_i = zeta^(i(gr+1))`,
    `gr = diag(0,..,4)`. Clashes checked against NOTATION tables 1-8: `calY` vs the uniformised `Y` of s08 (different section, other
    font); `T_i` vs sans `T_k` (QRR power sums) and `calT_k` (model): prefer writing `zeta^(i(gr+1))`; `lambda_m` vs `lambda_i` (torus
    weights) and the local `lambda` of s09: keep it inside one proof.
11. **Superseded wording.** UR's "unique up to `P -> P c(z)`" and GI s.1's "it can be taken graded" are now theorems (5.1 (b), 6.2);
    DB s.7 G1-G9 are closed except the dependency of item 5 on R1 (point 2).

## 12. Dependencies of each result (labels of `LABELS.md` and the new ones above)

* `lem:gau-comm`: none ((H) only for part (f)). `lem:gau-lift`: `lem:gau-comm` (b)-(f). `lem:gau-scal`: `lem:gau-comm` (c)-(f); (H), (EP).
* `thm:gau-U` (a): reduction 5.2, `lem:gau-lift`, `lem:gau-scal` (b), (c). (b): `lem:gau-lead` (from `lem:gau-comm` (c), (f) and
  `lem:gau-scal` (c)). (c): `lem:gau-degree`, (b). (d): the proofs of (a)-(c) (rational operations only). Hypotheses for the
  application: `prop:gau-quintic`, `prop:gau-equalprod` (outside this cluster).
* `prop:gau-graded`: `lem:gau-degree`, `thm:gau-U` (a), (b).
* `def:gau-can`, colour covariance of `P_can`: `prop:gau-graded`, `thm:gau-U` (b).
* `lem:gau-diag` step 4: GI s.5 steps 1-3 [W] or `lem:gau-lead`; (N0) from `def:set-frames`, `def:set-tau`, `prop:gau-quintic`,
  `prop:gau-equalprod`.
* `thm:gau-G` Cor. (a): (G) (`thm:gau-G`, i.e. `lem:gau-R1`, `lem:gau-vanmod`, `lem:gau-vanshape`, `lem:gau-vanq`, `lem:gau-uniq`),
  `thm:gau-U` (b), `def:gau-can`; `lem:gau-vanq` (c) for `Y^C_Q(0) = -e_1`.
* `cor:gau-colours`: `thm:gau-R0` (iv) (definition of `D`), `lem:set-canon` (exponents), the colour covariance of `P_can`,
  `lem:gau-colmod` (`prop:set-modrows` (b) uniqueness, `def:set-M0`), `lem:gau-colrows` (`prop:set-cyclicLR` shape [open, G1],
  `thm:set-mirror` (M) [G1], `thm:set-qrr` colour-`i` form [(QRR); draftnote (a)], `lem:gau-wkb` over `Q(zeta)`, `lem:set-hrot`,
  `D u_i = zeta^i L` (`def:set-genus0`), `def:set-rows` conventions `Delta_i^(1/2) = zeta^(4i) C_0/L`, `lem:gau-confexp`).

## 13. Files in this folder

`GAPS_G2.md` .. `GAPS_G2_6.md` (this note); `g2_core.py` (series class, PF4 Frobenius basis and log obstruction, the two systems);
`g2_linalg.py`, `g2_lemma2.py`, `g2_gauge.py`, `g2_G.py`, `g2_unique.py`, `g2_colours.py` and their logs `*.log` (s.10, with the
arguments listed there). Total running time a few minutes on one process (`g2_unique.py` is the longest, about one minute).
