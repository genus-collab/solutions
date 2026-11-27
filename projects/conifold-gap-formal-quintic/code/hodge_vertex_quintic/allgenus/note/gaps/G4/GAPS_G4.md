# Gap cluster G4 (section s06, the model): resolutions (key gaps/G4, 2026-11-15)

Author-facing working note of the gap-closing phase; not part of the note; nothing here enters a `.tex` file as is. AI-written and
AI-checked only; no human has read it. Labels: PROVED (written argument here, or a quoted theorem whose statement was read in the
full text, with hypotheses checked), VERIFIED (exact computation; range and field stated), CONJECTURAL, SKETCH, OPEN.
Pieces: this file (s.0 results, s.1 conventions, s.2 files); `GAPS_G4_2.md` (item 1, `prop:mod-Rtr`); `GAPS_G4_3.md` (item 2,
`lem:mod-xfun`); `GAPS_G4_4.md` (item 3, `thm:mod-Lreg`, `lem:mod-abdks`, `prop:mod-period`, `lem:mod-vanper`); `GAPS_G4_5.md`
(item 4, base `h = 1` of `thm:mod-V`; findings for the editor).
Abbreviations: DNOPS = [Dunin-Barkowski:2015caa] (arXiv:1509.06954, current arXiv PDF); EO = [Eynard:2007kz] (arXiv:math-ph/0702045v4,
7 May 2007, the PDF); ABDKS = [ABDKS25] (arXiv:2212.00320, current arXiv PDF); DOSS = arXiv:1211.4021. Project sources as in the design
(TRM = `allgenus/tr/TR_MODEL.md`, TRRF its referee, XY/XYR = `allgenus/xyswap/XYSWAP_MODEL.md` and referee, LL/LL2/LR = Stage-2 logleg
files, LG = `natural_form/lg/LG_MODEL.md`, UR = `natural_form/theory/UNIFORM_REDUCTION.md`).

## 0. Results

| item | label | resolution | in one line |
|---|---|---|---|
| 1 | `prop:mod-Rtr` (D6, G2) | **PROVED** | (2-31) of DNOPS, computed on `x = kappa W_v` by one integration by parts, IS the matrix of the five normalised model integrals `I_m`, `m = -2..2`: `X(zeta) = M_0(t)^T`, `t = -1/(kappa varsigma^4 zeta)`, exactly and to all orders; with the symplecticity of `M_0` (proved here, one paragraph) `R_TR(z) = M_0(z/(kappa varsigma^4))`; (2-32) holds identically and gives `M_0^-1 1hat`. Transfer `kappa = 4 -> 1` written out. No uniqueness theorem for homogeneous R-matrices and no sign-of-`kappa` argument is needed. |
| 2 | `lem:mod-xfun` (G3) | **PROVED** | Self-contained: the trace of `omega_(h,n)` over the five sheets of `x` vanishes (computed from the recursion: the trace of the kernel is `dX [1/(X - x(q)) - 1/(X - x(sigma q))] = 0`); then `Res_(varsigma_b) f(x) omega = 0` for EACH `b`, by the trace formula for residues and the distinct critical values. Also pinned in the literature: EO v4 Thm 4.4, (4-23), and the proof of Cor. 4.1, eq. (A-25), which states exactly this deduction. Precision: the pairing must be the residue pairing; the circle pairing of a leg with a pole at `p = 0` differs. |
| 3 | `thm:mod-Lreg` (+ `lem:mod-abdks`, `prop:mod-period`, `lem:mod-vanper`) (W5) | **PROVED** (mod ABDKS Thm 7.1, already in (Q); statement read in the full text) | Every step written: structure `omega = N prod dp_i/prod (p_i^5 - v)^K`, `N` in `Q[v][p]` (divisibility in `Q[v][p]`, since `p_i`, `p_i - p_j` are monic); pairing `= (-1)^n Res_oo^n`; expansion at `oo`. W5 dissolves: Thm 7.1 applies to `x = W_v` itself (any rational `x`); the transfer to XY's `x = z^4 + V/z` is one scaling, written out. Route 2 of the Gap Theorem: the "Borel identity" is replaced by an exact finite identity (`S_g = 4^g A_g`, the binomial sum and Lagrange inversion of route 1), so no Borel/Watson argument is needed. |
| 4 | `thm:mod-V` base `h = 1` (D5, G1) | **PROVED** (finite exact derivation by hand; no quoted theorem) | A local formula for `Int_r omega_(1,1)` (r a ray) from the definition of the recursion; at `varsigma_0 = 1`, `v = 1`: contribution `-1/60` per branch point, five equal contributions (rotation), total `-1/12`; homogeneity gives `<-log p>_1 = -1/(12 v)`. Two independent cross-checks written out: the vanishing period at `g = 1` (route of item 3), and the closed form of `omega_(1,1)` integrated by `u = p^5`. |

Evidence (all exact, all with failing controls, logs in this folder): `g4_rtr.py` (item 1, `Q(zeta5)`, through `zeta^-10`, `kappa =
1, 4, -1`, `varsigma = 1, 2`); `g4_w11.py` (items 2, 4); `g4_period.py` (item 3, `g <= 5`, `S_g` for `g <= 15`). Evidence only:
every resolution above is a written argument.

## 1. Conventions used in G4_2..G4_5 (those of NOTATION, NOTATION_2; local letters declared where used)

* Curve `x = W_v(p) = p^4/4 + v/p`, `y = p`, `omega_(0,2) = B = dp_1 dp_2/(p_1 - p_2)^2`, `v != 0`; branch points `varsigma_b =
  zeta^-b varsigma` (`zeta = e^(2 pi i/5)`, `varsigma^5 = v`), `x''(varsigma_b) = 5 varsigma_b^2`, critical values `x(varsigma_b) = (5/4)
  zeta^b varsigma^4` pairwise distinct. For the transfer: `x_kappa := kappa W_v` (`kappa = 4`: DNOPS's `f = 4 W_v`, `a = (0,0,0,0,4v)`).
* Kernel: engine convention `conv:mod-kernel`, `K(p_0, q) = (1/2) Int_(sigma q)^q B(p_0, .)/((y(q) - y(sigma q)) dx(q))`. The printed
  kernels of EO (4-3), (4-9) (`(1/2) Int_q^(qbar) B/((y(q) - y(qbar)) dx(q))`) and of ABDKS (1.1) with `omega_(1,0) = -y dx` are EQUAL to
  each other and opposite to the engine's; hence `omega_(h,n) = (-1)^n W^EO_(h,n) = (-1)^n omega^ABDKS_(h,n)` (checked on the PDFs).
* `sigma_b`: the local involution of `x` at `varsigma_b`; `frak r`: a ray from 0 to oo avoiding the branch points.
* Pairings: residue pairing `<phi, w> := Sum_b Res_(varsigma_b) phi w`; log leg `<-log p, w> := Int_(frak r) w`.
* Model data: `nu = (0,-1,3,2,1)`, `eps = (1,-1,-1,1,1)`, `Fs_b(t)` (T-system solution, `prop:set-modrows`), `Rhat^mod_b^(i)(t) =
  zeta^(i(b+1)) Fs_b(zeta^-i t)`, `Psi_0 = Rhat^mod(0)`, `M_0(t) = Psi_0^-1 Rhat^mod(t)`, `eta(H^a, H^b) = 5 delta_(a+b = 3)`.
  Note `(Psi_0)_(bi) = eps_b zeta^(i(b+1))` and `(Psi_0^-1)_(ib) = (1/5) eps_b zeta^(-i(b+1))` (since `eps_b^2 = 1` and
  `Sum_b zeta^((j-i)(b+1)) = 5 delta_ij`).
* DNOPS objects: local coordinates `w_j` with `x_kappa = x_kappa(varsigma_j) - w_j^2/2` (DNOPS s.2.6), `Delta_j` (2-30), `R^-1(zeta^-1)^j_i`
  (2-31), (2-32), `xi^i` (2-34), (2-35); the R-matrix action of DNOPS s.2.3 (leaves `R^-1(psi) e`, dilaton leaves `-psi (I - R^-1(psi)) e_1`,
  edges `(I (x) I - R^-1(psi') (x) R^-1(psi''))/(psi' + psi'')` applied to `eta^(ab) e_a (x) e_b`): this is the convention of
  [PPZ15] Def. 2.2; DNOPS Remark 2.1: Givental's is its transpose. Below `Lambda_kappa := -5 kappa` (a NUMBER, local to G4_2; not the
  leg map); a fixed square root `sqrt(Lambda_kappa)` is used for all `j` (uniform orientation, see G4_2 s.1).
* The Laplace variable of DNOPS is `zeta` in (2-31); to avoid the clash with the root of unity it is written `lam` in G4_2.

## 2. Files in this folder

| file | content |
|---|---|
| `GAPS_G4.md`, `_2` .. `_5` | this note |
| `g4_field.py` | exact `Q(zeta5)` arithmetic, power series (sqrt, inverse, reversion), Bernoulli numbers |
| `g4_rtr.py` -> `g4_rtr_10.log` | item 1: T-system `M_0`; (2-31) from local expansions with the finite part; `X = M_0(t)^T`; IBP form; (2-32); symplecticity; controls |
| `g4_w11.py` -> `g4_w11.log` | items 2, 4: local formula for `omega_(1,1)` and `<-log p>_1`; local formula = closed form; `Res f(x) omega = 0`; controls |
| `g4_period.py` -> `g4_period_5.log` | item 3: `S_g = -4^g B_2g/(2g)` (g <= 15); `Phi_g(0)`, `Phi_g(oo)` from ABDKS (7.10), `g <= 5`, `V = 1, 3`; controls |
| `pdftext.py` | reading aid (stdlib PDF text extraction) used to read the three papers; not evidence |

Run: `py -I -B g4_rtr.py 10`, `py -I -B g4_w11.py`, `py -I -B g4_period.py 5` (working directory = this folder or the scratch copy);
run times 10 s, 1.4 s, 0.2 s; one process each; logs end with "ALL CHECKS PASSED".
