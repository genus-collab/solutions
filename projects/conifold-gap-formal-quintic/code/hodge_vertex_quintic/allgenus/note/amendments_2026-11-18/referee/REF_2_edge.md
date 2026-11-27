# REF_2 -- referee report on `prop:mod-Rtr` (iii), `eq:mod-edge` (the edge condition)

Scope: `sections/s06_model_2.tex` l.52-57 (statement), `sections/s06_model_3.tex` l.21-25 (proof, Steps A-C),
l.29, l.34 (proof of `thm:mod-tr`, Step 1), `sections/s01_intro_3.tex` l.9 (Q), l.15 (readings).
Read-only on the note; nothing edited. Scripts (standard-library Python, double precision) in
`<user home>\AppData\Local\Temp\claude\d7c1b608d03d194a\37f24460-3601-4c85-b123-5156e19eb475\scratchpad\ref2\`
(`edge_check.py`, `mag.py`, `w11.py`).

## Verdict: SOUND AFTER FIXES

The identity `eq:mod-edge` is true, the proof (Steps A-C) is complete and correct, it is the right hypothesis in the
normalisation in which the note applies the quoted theorem, and `thm:mod-tr` Step 1 and (Q) use exactly what is
proved. No mathematical error and no gap in the proof. The fixes are three wording corrections (one of them a
misdescription of the misprint, repeated in three places) and one missing normalisation pin.

## 1. What was checked, and how

### 1.1 By computation (independent of the audit's scripts; no wolframscript)
`edge_check.py`: for (kappa, v) = (1,1), (4,1), (4, 0.7+0.4i), (1, -2.3+1.1i), (2.5-i, 0.3-0.8i): the local
coordinates w_j at the five critical points (x_kappa = x_kappa(c_j) - w^2/2, dw_j/dp(c_j) = sqrt(-5 kappa) c_j, one
root) by series reversion to order 12; R_TR^{-1} from its definition (-psi L_j[d xi^i], G_{-2} = -theta at i = j);
B^reg_ij as d_w d_w' log(p_i(w) - p_j(w')) (for i = j: log of the divided difference), as a bivariate series.
* `eq:mod-edge`: all 25 pairs (i,j), all coefficients psi_1^m psi_2^n with m + n <= 5: max residual 2e-16 ... 1e-15
  (the coefficients involved are of size 1e-3 ... 0.5, `mag.py`, so the test is not vacuous).
* R_TR^{-1}(0) = 1 (residual 1e-16), zero residue of d xi^i at its own branch point, (2-32) in the form of part (ii)
  (residual <= 1.4e-15), symplecticity R^{-1}(psi) R^{-1}(-psi)^T = 1 through degree 5 (2e-15).
* Step A as an identity of rational functions at 20 random points (relative residual 2e-15).
* The printed DNOPS (2-33) fails: for i != j the two sides differ by the SIGN (sample: coefficient of psi_1 at
  (i,j) = (0,1), kappa = 1, v = 1: left side -0.11708+0.08507i = minus the printed right side); max residual 1.1.
* `w11.py` (normalisation of the whole framework): omega_{1,1} computed from the recursion of `conv:mod-kernel`
  (residues at the five branch points) agrees with the graph formula (2-35) with: vertex Delta_i^{1/2}/24, leg
  R^{-1}_1, dilaton leaf of PPZ sign, and the loop with edge coefficient -(R^{-1}_1)^i_i = a_00 of B^reg_ii and
  the factor 1/2: ratio 1 to 1e-15 at four parameter sets (kappa = 1 and 4). This confirms that in the coordinates
  w_i and with the note's kernel the edge enters with the normalisation of `eq:mod-edge` (no factor 2, sign +), at
  the order psi^0 psi^0 of the diagonal.

### 1.2 By hand (re-derived)
* Operator D and the displayed bracket: d_1(f(p_1)/(p_1-p_2)^2) = f'(p_1)/(p_1-p_2)^2 - 2f(p_1)/(p_1-p_2)^3,
  d_2(f(p_2)/(p_1-p_2)^2) = f'(p_2)/(p_1-p_2)^2 + 2f(p_2)/(p_1-p_2)^3. Correct.
* Diagonal limit: with p_1 = p_2 + e the numerators are 2f' + f''e + f'''e^2/2 and -2f'e - f''e^2 - f'''e^3/3;
  the bracket tends to f'''(1/2 - 1/3) = f'''/6. Correct.
* At p_1 = c_k: the simple poles 2r_k/((c_k-p_2)^3 (p_1-c_k)) and -2r_k/((p_1-c_k)(c_k-p_2)^3) cancel; principal
  part -r_k/((p_1-c_k)^2 (c_k-p_2)^2); (dp/dw_k)^2 = -1/x''(c_k) = -r_k. Correct. Holomorphy at 0 and infinity
  (the bracket is O(p_1^{-2}) at infinity): correct. Uniqueness on P^1: correct.
* Step B: F = Omega/dx_1 has the single polar term w^{-1}; `lem:mod-laplace` (i) applies (F dx = Omega has no
  residue), dF has only the polar term w^{-2} (no w^{-1}), so L_i (x) L_j [D Omega] is defined. Correct.
* Step C, i = j: (w^{-2} + w'^{-2})/(w-w')^2 - 2/(w w'(w-w')^2) = 1/(w^2 w'^2); L (x) L of dw dw'/(w^2 w'^2) is
  G_{-2}(theta_1) G_{-2}(theta_2) = theta_1 theta_2; the double Laurent expansion of D B^reg_ii with poles only on
  the coordinate axes is unique, so the termwise transform is legitimate. Signs and the division by
  -theta_1 theta_2: correct. No factor sqrt(2 pi) enters: the statement and proof use only the formal L.
* Index convention: R^{-1}(psi)^j_k = -psi L_j[d xi^k], upper index = Laplace point = row; the right side of
  `eq:mod-edge` is (1 - R^{-1}(psi_1) R^{-1}(psi_2)^T)^{ij}, i.e. (psi_1+psi_2) times the PPZ edge in the
  orthonormal basis. DOSS (4-1) has sum_s R^i_s(-z) R^j_s(-w) with the upper index the column (DOSS Remark 2.1)
  and R(-z)^T = R(z)^{-1}: the same matrix. Consistent with parts (i), (ii) and with `sec:sgn-R`.
* Normalisation (the question "is this the hypothesis the theorem needs"): in the coordinates w_i the half-edge
  forms have leading term (2d+1)!! w^{-2d-2} dw (d((d/dx)^d xi^i), xi^i = -1/w_i + ...), which are the forms of
  DOSS (3-55); the residue pairing (3-56) of B^reg with two of them is a_{2d_1,2d_2} (2d_1-1)!! (2d_2-1)!!, i.e.
  the coefficient of psi_1^{d_1} psi_2^{d_2} of L (x) L [B^reg]. So `eq:mod-edge` is exactly "TR edge weight = PPZ
  edge of R_TR". The edge weight is insensitive to a common sign of the half-edge forms (two half-edges per edge).

### 1.3 By reading only
DOSS (3-1), (3-2), (3-55)-(3-57), (3-60)-(3-66), Thm 4.1 with (4-1)-(4-3) and DNOPS 2.6 (page image p.12) in the
stored extracts under `audit_2026-11-18/work/L3/`. The agreement of the edge normalisation at orders beyond
psi^0 psi^0 with the quoted theorem was checked only by the reading in 1.2 (last item), not by a higher-genus
computation (that would need omega_{1,2} or omega_{0,5}).

## 2. Findings

### F1 (wording; three places) -- the misprint in DNOPS (2-33) is misdescribed
The printed (2-33) reads L (x) L [B] = + sum_k R^{-1}(psi_1)^i_k R^{-1}(psi_2)^j_k / (psi_1 + psi_2). The correct
statement is (delta_ij - sum_k ...)/(psi_1+psi_2) for B^reg. So the print differs by the term delta_ij AND by the
sign of the sum; for i != j, where no delta occurs, it is still false (wrong sign; confirmed numerically, 1.1).
(If the double integral over the full B is regularised, the correct right side is - sum_k ... /(psi_1+psi_2): again
the sign.) "Omits the term delta_ij" would leave a reader to restore delta_ij + sum, which is wrong.

* `sections/s06_model_2.tex` l.57.
  Old: `The form printed in \cite[(2-33)]{Dunin-Barkowski:2015caa} omits the term $\delta_{ij}$ and is not used.`
  New: `The form printed in \cite[(2-33)]{Dunin-Barkowski:2015caa} has $+\sum_k$ in place of $\delta_{ij}-\sum_k$ (so for $i\neq j$ it differs from \eqref{eq:mod-edge} by a sign) and is not used.`
* `sections/s06_model_3.tex` l.29.
  Old: `(the printed \cite[(2-33)]{Dunin-Barkowski:2015caa} omits the term $\delta_{ij}$)`
  New: `(the printed \cite[(2-33)]{Dunin-Barkowski:2015caa} has $+\sum_k$ in place of $\delta_{ij}-\sum_k$)`
* `sections/s01_intro_3.tex` l.9.
  Old: `is misprinted (it omits the term $\delta_{ij}$) and is not used`
  New: `is misprinted (it has $+\sum_k$ in place of $\delta_{ij}-\sum_k$) and is not used`

### F2 (wording / pin) -- the coordinate in which DOSS (4-1), (3-57) is read is not stated
DOSS (3-1) normalises the local coordinate by x = z^2 + a_i; the note uses the coordinates of DNOPS 2.6,
x = x(c_i) - w_i^2/2. The weights (3-55), (3-57) (double factorials without powers of 2) are those of the
coordinate w_i (see 1.2), and DOSS's own Section 3 is not uniform in this respect ((3-61) with x = z^2 would give
(2d-1)!!/2^d). The sentence "This is the edge condition [(4-1), with (3-2) and (3-57)]" is correct only with the
coordinate pinned; as written a reader cannot tell which normalisation is meant. This is the one place where the
note interprets, rather than quotes, the hypothesis of (Q); it should say so in one sentence.

* `sections/s06_model_2.tex` l.57.
  Old: `This is the edge condition \cite[(4-1), with (3-2) and (3-57)]{DuninBarkowski:2012bw}: the expansion of $\omega_{0,2}$ at the branch points is the edge of $R_{\mathrm{TR}}$.`
  New: `This is the edge condition \cite[(4-1), with (3-2) and (3-57)]{DuninBarkowski:2012bw}, written in the coordinates $w_i$ of \cite[\S2.6]{Dunin-Barkowski:2015caa} (the coordinate of \cite[(3-1)]{DuninBarkowski:2012bw} is normalised by $x=z^2+a_i$ instead): the coefficient of $\psi_1^{d_1}\psi_2^{d_2}$ in $(\mathsf L_i\otimes\mathsf L_j)[B^{\mathrm{reg}}_{ij}]$ is $(2d_1-1)!!\,(2d_2-1)!!\,a_{2d_1,2d_2}$, the edge weight (3-57), and the right side divided by $\psi_1+\psi_2$ is the edge of \cite[Def.~2.2]{PPZ15} for $R_{\mathrm{TR}}$ in the orthonormal basis $e^{\mathrm{TR}}_j$ (\cref{sec:sgn-R}). So the expansion of $\omega_{0,2}$ at the branch points is the edge of $R_{\mathrm{TR}}$.`
* `sections/s01_intro_3.tex` l.15 (list of readings), end of the paragraph.
  Old: `and the value $\mathsf G_{-2}=-\vartheta$ for the polar term $w^{-2}\dd w$ (cf.\ \cite[(3-64)]{DuninBarkowski:2012bw}).`
  New: `and the value $\mathsf G_{-2}=-\vartheta$ for the polar term $w^{-2}\dd w$ (cf.\ \cite[(3-64)]{DuninBarkowski:2012bw}); and the edge condition \cite[(4-1)]{DuninBarkowski:2012bw} read in the coordinates $w_i$ of \cite[\S2.6]{Dunin-Barkowski:2015caa}, that is, as \eqref{eq:mod-edge} (\cref{prop:mod-Rtr}~(iii)).`

### F3 (wording, optional) -- "symmetric bidifferential" in the definition of D
`sections/s06_model_3.tex` l.21. D is applied to the expansion of omega_{0,2} at (c_i, c_j), i != j, and to
d xi^k d xi^k in two different local coordinates; symmetry is neither available there nor used.
  Old: `for a symmetric bidifferential $\Omega$, $\mathsf D\Omega:=`
  New: `for a bidifferential $\Omega$, $\mathsf D\Omega:=`

### F4 (wording, optional) -- one clause in Step C to make the termwise transform visibly legitimate
`sections/s06_model_3.tex` l.25.
  Old: `So $\mathsf DB^{\mathrm{reg}}_{ii}=\sum_k\dd\xi^k\,\dd\xi^k-\dd w\,\dd w'/(w^2w'^2)$, and by Step~B`
  New: `So $\mathsf DB^{\mathrm{reg}}_{ii}=\sum_k\dd\xi^k\,\dd\xi^k-\dd w\,\dd w'/(w^2w'^2)$, an identity of double Laurent series in $w,w'$ (both sides have poles only on $w=0$ and $w'=0$, and $(\mathsf L_i\otimes\mathsf L_i)$ of a product $\alpha(w)\beta(w')$ is $\mathsf L_i[\alpha]\,\mathsf L_i[\beta]$), and by Step~B`

## 3. Consequences elsewhere (checked by reading and grep over `sections/*.tex`, `main.tex`)
* `thm:mod-tr`, Step 1 (`s06_model_3.tex` l.29, l.34): the three hypotheses are stated as symplecticity, (2-32),
  and `eq:mod-edge`, and discharged by `prop:mod-Rtr` (i), (ii), (iii) at kappa = 4; (iii) is proved for every
  kappa != 0, so this is exactly what is proved. (Symplecticity is now also a consequence of (iii); harmless.)
* (Q), `s01_intro_3.tex` l.9: lists DOSS Thm 4.1 with (4-1)-(4-3), the weights (2-22)-(2-26), DNOPS 2.6
  (2-30)-(2-32), (2-34), (2-35) -- (2-33) excluded -- and says the hypotheses are proved in `prop:mod-Rtr`.
  Correct after F1.
* No remaining reliance on DNOPS (2-33), Prop. 4.3, the sentence before Thm 8.1, or (8-1): the only occurrences
  are (a) the three "is not used" remarks; (b) `s06_model_3.tex` l.34, (8-1) "for orientation" and Prop. 4.3 "is
  not used" (the remark that it is equivalent to (2-32) agrees with DNOPS Thm 4.1 / Cor. 4.2, read in the stored
  extract); (c) the table rows `b_signs.tex` l.52, l.58 ((8-1) as a dictionary entry, "not used"). None is
  load-bearing.
* `b_signs_2.tex` l.25 (table row "Laplace transforms") mentions only (2-31), (2-32); no change needed, but the row
  could add "and the edge condition `eq:mod-edge`" for completeness (optional).

## 4. Summary table
| # | severity | location | status |
|---|----------|----------|--------|
| F1 | wording (misdescription of a source, 3 places) | s06_model_2 l.57; s06_model_3 l.29; s01_intro_3 l.9 | fix given |
| F2 | wording / normalisation pin | s06_model_2 l.57; s01_intro_3 l.15 | fix given |
| F3 | wording, optional | s06_model_3 l.21 | fix given |
| F4 | wording, optional | s06_model_3 l.25 | fix given |
No finding of severity "error" or "gap" in the mathematics of (iii) or its proof.
